//! Initializes and runs a simulation

use crate::boundary_conditions::*;
use crate::constants::Constants;
use crate::lbm_algorithm::*;
use crate::simulation_parsing::*;
use crate::vessels::*;
use crate::vessels_parsing::*;
use itertools_num::linspace;
use nalgebra::{Matrix2, Vector2};
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
    sm_soluce: SlotMap<VesselKey, f64>,
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
        let mut sm_soluce = SlotMap::with_key();

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
                vessels[i].rhs_f1 = sm_soluce.insert(0.0);
                vessels[i].lhs_f2_give = sm_soluce.insert(0.0);
                vessels[j].lhs_f2_recv = vessels[i].lhs_f2_give;

                let cp2 =
                    1.5 * vessels[i].cells.gamma[vessels[i].x_last] * vessels[i].cells.A0[vessels[i].x_last].sqrt();
                let cp2_r = 1.5 * vessels[j].cells.gamma[0] * vessels[j].cells.A0[0].sqrt();
                let coef = vessels[j].cells.A0[0] * cp2 / (vessels[i].cells.A0[vessels[i].x_last] * cp2_r);
                vessels[i].rhs_coef = ((vessels[i].cells.beta[vessels[i].x_last]) / (vessels[j].cells.beta[0])).powi(2)
                    * ((vessels[j].cells.A0[0]) / (vessels[i].cells.A0[vessels[i].x_last]));
                // debug!("{} {}", vessels[i].rhs_coef, coef);
                // vessels[i].rhs_coef = coef;

                let lhs_child_key = vessels[j].lhs_give;
                let rhs_parent_key = vessels[i].rhs_give;

                vessels[i].rhs_recv.push(lhs_child_key);
                vessels[j].lhs_recv.push(rhs_parent_key);
            }

            if vessels[i].children.len() == 0 {
                vessels[i].rhs_f1 = sm_soluce.insert(0.0);
            }

            if vessels[i].is_inlet {
                vessels[i].lhs_f2_recv = sm_soluce.insert(0.0);
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
            sm,
            sm_soluce,
            inlet: Inlet::new(&inlet_data),
            file: File::create("res/A").unwrap(),
            save_nb: 0,
        }
    }

    pub fn run<T: Compute>(&mut self) {
        let compute = T::new(self.consts.clone());

        for mut v in self.vessels.iter_mut() {
            let A_lhs = match v.parent_nb {
                1 => None, //Some(self.sm[v.lhs_recv[0]].A),
                _ => None,
            };

            let A_rhs = match v.children.len() {
                1 => None, //Some(self.sm[v.rhs_recv[0]].A),
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
                    info!("{}: {:?}", v.name, &v.cells.A[..5]);
                    // info!("A: {:?}", &v.cells.A[100..105]);
                    info!("{}: {:?}", v.name, &v.cells.A[v.x_last - 5..v.x_last]);
                    let mut A_file = &v.A_file;
                    let mut u_file = &v.u_file;

                    // let A: Vec<f64> = v.cells.A.iter().step_by(1).map(|&x| x).collect();
                    let A: Vec<f64> = v.cells.A.iter().step_by(1).map(|&x| x).collect();
                    // let A: Vec<f64> = zip_map!(v.cells, |(A, u)| A + u); //v.cells.A.iter().step_by(1).map(|&x| x).collect();
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
                1 => Some(v.cells.A[0]), //Some(self.sm[v.lhs_recv[0]].A),
                _ => None,
            };

            let A_rhs = match v.children.len() {
                1 => Some(v.cells.A[v.x_last]), //Some(self.sm[v.rhs_recv[0]].A),
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
            match v.outflow {
                Some(Outflow::NonReflective) => (),
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
                    let f2 = v.cells.f2[v.x_last - 1]; // Assume already streamed
                    let F = ((self.consts.dt * aL) / (2.0 * self.consts.c));
                    let f1 = F + f2 - (Qn / self.consts.c);

                    self.sm_soluce[v.rhs_f1] = f1;
                }

                None => (),
            }
        }

        for v in self.vessels.iter() {
            match v.is_inlet {
                true => {
                    let mut u = self.inlet.flow_at_time(self.current_time) / v.cells.A[0];
                    if self.current_time > 1.1 {
                        u = 0.0;
                    }

                    let A = (2.0 * v.cells.f1[1] + v.cells.f0[0]
                        - (v.cells.F[0] * self.consts.dt) / (2.0 * self.consts.c))
                        / (1.0 - (u / self.consts.c));

                    self.sm_soluce[v.lhs_f2_recv] = v.cells.f1[1] + A * (u / self.consts.c)
                        - (v.cells.F[0] * self.consts.dt) / (2.0 * self.consts.c);
                }

                false => (),
            }

            match v.outflow {
                Some(_) => (),
                None => match v.children.len() {
                    1 => {
                        let rhs = self.get_rhs_boundary(v);
                        let b1 = rhs[0].f.0 + rhs[0].f.1
                            - v.rhs_coef * (v.cells.f0[v.x_last] + v.cells.f2[v.x_last] - v.cells.A0[0])
                            - rhs[0].A0;
                        let b2 = v.cells.f2[v.x_last]
                            + rhs[0].f.1
                            + (self.consts.dt / (2.0 * self.consts.c)) * (v.cells.F[v.x_last] - rhs[0].F);

                        let A = Matrix2::new(v.rhs_coef, -1.0, 1.0, 1.0);
                        let b = Vector2::new(b1, b2);
                        let x = A.lu().solve(&b).unwrap();

                        self.sm_soluce[v.rhs_f1] = x[0];
                        self.sm_soluce[v.lhs_f2_give] = x[1];
                    }
                    _ => panic!("Vessel {} has {} children", v.name, v.children.len()),
                },
            }
        }

        for mut v in self.vessels.iter_mut() {
            v.cells.f2.pop_back().unwrap();
            v.cells.f2.push_front(self.sm_soluce[v.lhs_f2_recv]);

            v.cells.f1.pop_front().unwrap();
            v.cells.f1.push_back(self.sm_soluce[v.rhs_f1]);
        }

        self.current_time += self.consts.dt;
        self.time_since_last_save += self.consts.dt;
        self.current_iter += 1;
    }

    pub fn get_rhs_boundary(&self, v: &Vessel) -> Vec<VesselBoundary> {
        v.rhs_recv.iter().map(|key| self.sm[*key]).collect()
    }
}
