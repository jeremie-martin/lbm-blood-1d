//! Initializes and runs a simulation

use crate::boundary_conditions::*;
use crate::constants::Constants;
use crate::lbm_algorithm::*;
use crate::simulation_parsing::*;
use crate::vessels::*;
use crate::vessels_parsing::*;
use itertools_num::linspace;
use serde::{Deserialize, Serialize};
use slotmap::{new_key_type, Key, SlotMap, Slottable};
use std::collections::VecDeque;
use std::fs::File;
use std::io::prelude::*;
use std::io::Write;
use tracing::{event, info, instrument, span, warn, Level}; // This trait adds methods to writeable types

#[derive(Debug, Copy, Clone)]
pub struct VesselBoundary {
    pub A0: f64,
    pub A: f64,
    pub F: f64,
    pub f: (f64, f64, f64),
}

/// Contains the simulation parameters and the vessels.
#[derive(Debug)]
pub struct Simulation {
    /// Implementation of the LBM steps
    // pub compute: Box<T>,
    pub current_time: f64,
    pub current_iter: u64,
    /// Duration of the simulation
    pub total_time: f64,
    /// Cardiovascular network
    pub vessels: Vec<Vessel>,
    /// Constants related to the cardiovascular network
    pub consts: Constants,
    /// Inlet (heart)
    pub inlet: Inlet,
    pub file: File,
    /// Save every `x` second
    pub time_between_save: f64,
    /// Elapsed time since last save
    pub time_since_last_save: f64,
    pub save_nb: u64,
    sm: SlotMap<VesselKey, VesselBoundary>,
}

impl Simulation {
    /// Returns a structure ready to be simulated, given the `path` to a .json describing both the simulation parameters and the cardiovascular network.
    pub fn new(path: &str, dx: f64, time_between_save: f64) -> Simulation {
        let mut parse = SimulationParsing::read_json(&path.to_string());
        info!("Unmarshalled {}", path);

        let inlet_data = SimulationParsing::read_inlet(&parse.inlet_path);
        info!("Unmarshalled {}", parse.inlet_path);

        let total_time = parse.total_time;

        let consts = Constants::new(&parse, dx);

        info!("Precomputed constants");
        info!(
            ">dt = {:.4e}, c = {:.4e}, cs = {:.4e}, cs2 = {:.4e}",
            consts.dt, consts.c, consts.cs, consts.cs2
        );
        info!(
            ">nu = {:.4e}, tau = {:.4e}, omega = {:.4e}",
            consts.nu, consts.tau, consts.omega
        );

        let mut vessels: Vec<Vessel> = parse.vessels.iter().map(|v| Vessel::new(&v, consts.clone())).collect();

        let mut sm = SlotMap::with_key();

        for mut v in vessels.iter_mut() {
            v.lhs_give = sm.insert(VesselBoundary {
                A0: v.cells.A0[0],
                A: v.cells.A[0],
                F: v.cells.F[0],
                f: (v.cells.f0[0], v.cells.f1[0], v.cells.f2[0]),
            });

            v.rhs_give = sm.insert(VesselBoundary {
                A0: v.cells.A0[v.x_last],
                A: v.cells.A[v.x_last],
                F: v.cells.F[v.x_last],
                f: (v.cells.f0[v.x_last], v.cells.f1[v.x_last], v.cells.f2[v.x_last]),
            });
        }

        for i in 0..vessels.len() {
            for &j in vessels[i].children.clone().iter() {
                vessels[j].parent_nb += 1;

                let lhs_child_key = vessels[j].lhs_give;
                let rhs_parent_key = vessels[i].rhs_give;

                vessels[i].rhs_recv.push(lhs_child_key);
                vessels[j].lhs_recv.push(rhs_parent_key);
            }
        }

        Simulation {
            consts,
            current_time: 0.0,
            current_iter: 0,
            total_time,
            vessels,
            time_between_save,
            time_since_last_save: 0.0,
            sm: sm,
            inlet: Inlet::new(&inlet_data),
            file: File::create("res/A").unwrap(),
            save_nb: 0,
        }
    }

    pub fn run<T: Compute>(&mut self) {
        let compute = T::new(self.consts.clone());

        for mut v in self.vessels.iter_mut() {
            let A_lhs = match v.parent_nb {
                1 => Some(self.sm[v.lhs_recv[0]].A),
                _ => None,
            };

            let A_rhs = match v.children.len() {
                1 => Some(self.sm[v.rhs_recv[0]].A),
                _ => None,
            };

            compute.forcing_term(&mut v, A_lhs, A_rhs);
        }

        for v in &mut self.vessels {
            zip_for_each!(v.cells, |(A, mut u, F)| { *u = compute.velocity_init(A, u, F) });

            compute.dev(v);
            zip_for_each!(v.cells, |(A, u, F, deriv, mut f0, mut f1, mut f2)| {
                let feq = compute.populations_init(A, u, F, deriv);
                *f0 = feq.0;
                *f1 = feq.1;
                *f2 = feq.2;
            });

            zip_for_each!(v.cells, |(mut A, f0, f1, f2)| { *A = compute.area((f0, f1, f2)) });
            zip_for_each!(v.cells, |(A, u, F, beta, A0, f0, f1, f2, mut stress)| {
                *stress = compute.stress(A, u, F, beta, A0, (f0, f1, f2))
            });

            match v.outflow {
                Some(Outflow::WK3(ref mut wk3)) => {
                    let P = v.cells.beta[v.x_last] * ((v.cells.A[v.x_last] / v.cells.A0[v.x_last]).sqrt() - 1.0);
                    let Q = v.cells.A[v.x_last] * v.cells.u[v.x_last];
                    wk3.P_old = VecDeque::from(vec![P; 2]);
                    wk3.Q_old = VecDeque::from(vec![Q; 2]);
                }
                _ => (),
            }

            self.sm[v.lhs_give].A = v.cells.A[0];
            self.sm[v.lhs_give].F = v.cells.F[0];
            self.sm[v.lhs_give].f = (v.cells.f0[0], v.cells.f1[0], v.cells.f2[0]);

            self.sm[v.rhs_give].A = v.cells.A[v.x_last];
            self.sm[v.rhs_give].F = v.cells.F[v.x_last];
            self.sm[v.rhs_give].f = (v.cells.f0[v.x_last], v.cells.f1[v.x_last], v.cells.f2[v.x_last]);
        }

        info!("Initialized vessels");

        while self.current_time < 3.3 {
            if (self.current_iter == 0 || self.time_since_last_save >= self.time_between_save) {
                for v in self.vessels.iter() {
                    info!("Current time: {:.6}s (iter {})", self.current_time, self.current_iter);
                    info!("A: {:?}", &v.cells.A[..5]);
                    // info!("A: {:?}", &v.cells.A[100..105]);
                    info!("A: {:?}", &v.cells.A[v.x_last - 5..v.x_last]);
                    let mut A_file = &v.A_file;
                    let mut u_file = &v.u_file;

                    let A: Vec<f64> = v.cells.A.iter().step_by(1).map(|&x| x).collect();
                    info!("A: {:?}", A.len());

                    unsafe {
                        A_file.write(std::slice::from_raw_parts(A.as_ptr() as *const u8, A.len() * 8));
                    }
                }

                self.time_since_last_save = 0.0;
                self.save_nb += 1;
            }
            self.one_step(&compute);
        }

        info!("Simulation done (tot: {} save)", self.save_nb);
    }

    pub fn one_step<T: Compute>(&mut self, compute: &T) {
        // Save

        for mut v in self.vessels.iter_mut() {
            let A_lhs = match v.parent_nb {
                1 => Some(self.sm[v.lhs_recv[0]].A),
                _ => None,
            };

            let A_rhs = match v.children.len() {
                1 => Some(self.sm[v.rhs_recv[0]].A),
                _ => None,
            };

            compute.forcing_term(&mut v, A_lhs, A_rhs);
        }

        for mut v in self.vessels.iter_mut() {
            zip_for_each!(v.cells, |(mut A, f0, f1, f2)| { *A = compute.area((f0, f1, f2)) });
            zip_for_each!(v.cells, |(A, mut u, F, f0, f1, f2)| {
                *u = compute.velocity(A, F, (f0, f1, f2))
            });

            zip_for_each!(v.cells, |(A, u, F, mut f0, mut f1, mut f2)| {
                let bgk = compute.BGK(A, u, F, (f0, f1, f2));
                *f0 = bgk.0;
                *f1 = bgk.1;
                *f2 = bgk.2;
            });

            self.sm[v.lhs_give].A = v.cells.A[0];
            self.sm[v.lhs_give].F = v.cells.F[0];
            self.sm[v.lhs_give].f = (v.cells.f0[0], v.cells.f1[0], v.cells.f2[0]);

            self.sm[v.rhs_give].A = v.cells.A[v.x_last];
            self.sm[v.rhs_give].F = v.cells.F[v.x_last];
            self.sm[v.rhs_give].f = (v.cells.f0[v.x_last], v.cells.f1[v.x_last], v.cells.f2[v.x_last]);
        }

        for mut v in self.vessels.iter_mut() {
            v.cells.f2.pop_back().unwrap();

            let f2 = match v.is_inlet {
                true => {
                    let mut u = self.inlet.flow_at_time(self.current_time) / v.cells.A[0];
                    if self.current_time > 1.1 {
                        u = 0.0;
                    }

                    let A = (2.0 * v.cells.f1[1] + v.cells.f0[0]
                        - (v.cells.F[0] * self.consts.dt) / (2.0 * self.consts.c))
                        / (1.0 - (u / self.consts.c));

                    v.cells.f1[1] + A * (u / self.consts.c) - (v.cells.F[0] * self.consts.dt) / (2.0 * self.consts.c)
                }

                false => match (v.parent_nb) {
                    1 => self.sm[v.lhs_recv[0]].f.2,
                    _ => panic!("Vessel {} has {} parents", v.name, v.parent_nb),
                },
            };
            // // TODO: wrap this map in a macro
            // let a = vessels[0].rhs_recv.iter().map(|key| sm[*key].f).collect::<Vec<f64>>(); // [2.0, 3.0]
            // let b = vessels[1].lhs_recv.iter().map(|key| sm[*key].f).collect::<Vec<f64>>(); // [1.0]
            // let c = vessels[2].lhs_recv.iter().map(|key| sm[*key].f).collect::<Vec<f64>>(); // [1.0]

            v.cells.f2.push_front(f2);

            v.cells.f1.pop_front().unwrap();

            let f1 = match v.outflow {
                Some(Outflow::NonReflective) => 0.0,
                Some(Outflow::WK3(ref mut wk3)) => {
                    let (P, Pn, Q, Qn) = wk3.owo(
                        v.cells.A[v.x_last],
                        v.cells.u[v.x_last],
                        v.cells.A0[v.x_last],
                        v.cells.beta[v.x_last],
                        v.consts.dt,
                    );

                    let A = v.cells.A0[v.x_last] * ((Pn / v.cells.beta[v.x_last]) + 1.0).powi(2);
                    let u = Qn / A;

                    wk3.P_old.pop_back();
                    wk3.P_old.push_front(P);
                    wk3.Q_old.pop_back();
                    wk3.Q_old.push_front(Q);

                    let aL = 2.0 * v.cells.F[v.x_last - 1] - v.cells.F[v.x_last - 2];
                    let f2 = v.cells.f2[v.x_last];
                    let F = ((self.consts.dt * aL) / (2.0 * self.consts.c));
                    let f1 = F + f2 - (Qn / self.consts.c);

                    f1
                }

                None => match (v.children.len()) {
                    1 => self.sm[v.rhs_recv[0]].f.1,
                    _ => panic!("Vessel {} has {} children", v.name, v.children.len()),
                },
            };

            v.cells.f1.push_back(f1);

            // let back = match v.outflow {
            //     Some(Outflow::NonReflective) => compute.non_reflective(v),
            //     Some(Outflow::WK3(_)) => compute.non_reflective(v),
            //     None => panic!("Vessel {} ({}) outlet type is not set", v.name, v.id),
            // };
        }

        self.current_time += self.consts.dt;
        self.time_since_last_save += self.consts.dt;
        self.current_iter += 1;
    }
}
