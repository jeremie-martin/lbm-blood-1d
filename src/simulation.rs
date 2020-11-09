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
use tracing::{event, info, instrument, span, warn, Level};

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

        Simulation {
            consts,
            current_time: 0.0,
            current_iter: 0,
            total_time,
            vessels,
            time_between_save,
            time_since_last_save: 0.0,
            inlet: Inlet::new(&inlet_data),
            file: File::create("res/A").unwrap(),
        }
    }

    pub fn run<T: Compute>(&mut self) {
        let compute = T::new(self.consts.clone());

        for v in &mut self.vessels {
            if v.is_inlet == false {
                continue;
            }

            info!(
                "f0: {} {} {} {} {}",
                &v.cells.f0[0], &v.cells.f0[1], &v.cells.f0[2], &v.cells.f0[3], &v.cells.f0[4],
            );
            info!(
                "f1: {} {} {} {} {}",
                &v.cells.f1[0], &v.cells.f1[1], &v.cells.f1[2], &v.cells.f1[3], &v.cells.f1[4],
            );
            info!(
                "f2: {} {} {} {} {}",
                &v.cells.f2[0], &v.cells.f2[1], &v.cells.f2[2], &v.cells.f2[3], &v.cells.f2[4],
            );
            compute.forcing_term(v);

            zip_for_each!(v.cells, |(A, mut u, F)| { *u = compute.velocity_init(A, u, F) });

            compute.dev(v);
            zip_for_each!(v.cells, |(A, u, F, deriv, mut f0, mut f1, mut f2)| {
                // println!("deriv {}", deriv);
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
        }

        while self.current_time < 3.3 {
            self.one_step(&compute);
        }

        info!("Initialized vessels");
    }

    pub fn one_step<T: Compute>(&mut self, compute: &T) {
        // Save

        let mut file = self.file.try_clone().unwrap();
        for mut v in self.vessels.iter_mut() {
            if v.is_inlet == false {
                continue;
            }
            if self.time_since_last_save >= self.time_between_save {
                info!("Current time: {:.6}s (iter {})", self.current_time, self.current_iter);
                info!("A: {:?}", &v.cells.A[..5]);
                // info!("A: {:?}", &v.cells.A[100..105]);
                info!("A: {:?}", &v.cells.A[v.x_last - 5..v.x_last]);
                zip_for_each!(v.cells, |(A, u, F)| {
                    // for A in &v.cells.A {
                    write!(file, "{} ", A);
                });
                write!(self.file, "\n");

                self.time_since_last_save = 0.0;
            }

            // info!("New iter!");
            compute.forcing_term(&mut v);
            // info!("F: {:?}", &v.cells.u[..5]);

            zip_for_each!(v.cells, |(mut A, f0, f1, f2)| { *A = compute.area((f0, f1, f2)) });
            // info!("A: {:?}", &v.cells.A[..5]);
            zip_for_each!(v.cells, |(A, mut u, F, f0, f1, f2)| {
                *u = compute.velocity(A, F, (f0, f1, f2))
            });

            // info!("u: {:?}", &v.cells.u[..5]);

            zip_for_each!(v.cells, |(A, u, F, mut f0, mut f1, mut f2)| {
                let bgk = compute.BGK(A, u, F, (f0, f1, f2));
                *f0 = bgk.0;
                *f1 = bgk.1;
                *f2 = bgk.2;
            });

            let backkk = v.cells.f2.pop_back().unwrap();
            let mut u = self.inlet.flow_at_time(self.current_time) / v.cells.A[0];
            if self.current_time > 1.1 {
                u = 0.0;
            }

            let A = (2.0 * v.cells.f1[1] + v.cells.f0[0] - (v.cells.F[0] * self.consts.dt) / (2.0 * self.consts.c))
                / (1.0 - (u / self.consts.c));

            // v.cells.f2.push_front(v.cells.f1[1] + A * (u / self.consts.c));
            v.cells.f2.push_front(
                v.cells.f1[1] + A * (u / self.consts.c) - (v.cells.F[0] * self.consts.dt) / (2.0 * self.consts.c),
            );

            let mut front = v.cells.f1.pop_front().unwrap();
            let idx = v.x_last;
            // let f1 = v.cells.f1[idx - 1];
            let mut QQ = 0.0f64;
            match v.outflow {
                Some(Outflow::NonReflective) => panic!("owo"),
                Some(Outflow::WK3(ref mut wk3)) => {
                    let idx = v.x_last;
                    let (P, Pn, Q, Qn) = wk3.owo(
                        v.cells.A[idx],
                        v.cells.u[idx],
                        v.cells.A0[idx],
                        v.cells.beta[idx],
                        v.consts.dt,
                    );

                    let A = v.cells.A0[idx] * ((Pn / v.cells.beta[idx]) + 1.0).powi(2);
                    let u = Qn / A;
                    QQ = Qn;

                    wk3.P_old.pop_back();
                    wk3.P_old.push_front(P);
                    wk3.Q_old.pop_back();
                    wk3.Q_old.push_front(Q);

                    let A0 = v.cells.A0[v.x_last];
                    let aL = 2.0 * v.cells.F[v.x_last - 1] - v.cells.F[v.x_last - 2];
                    let f2 = v.cells.f2[v.x_last];
                    let F = ((self.consts.dt * aL) / (2.0 * self.consts.c));
                    let f1 = F + f2 - (Qn / self.consts.c);

                    v.cells.f1.push_back(f1);
                }
                None => panic!("Vessel {} ({}) outlet type is not set", v.name, v.id),
            };

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
