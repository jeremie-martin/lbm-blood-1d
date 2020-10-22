//! Initializes and runs a simulation

use crate::boundary_conditions::Inlet;
use crate::constants::Constants;
use crate::simulation_parsing::*;
use crate::vessels::*;
use itertools_num::linspace;
use replace_with::replace_with;
use serde::{Deserialize, Serialize};
use splines::{Interpolation, Key, Spline};
use std::fs::File;
use tracing::{event, info, instrument, span, warn, Level};

/// Contains the simulation parameters and the vessels.
#[derive(Debug)]
pub struct Simulation {
    /// [Constants](crate::Constants)
    pub consts: Constants,
    /// Duration of the simulation
    pub total_time: f64,
    /// [Vessels](crate::Vessel)
    pub vessels: Vec<Vessel>,
    /// Uniformly sampled inflow (w.r.t. dt)
    pub inlet: Inlet,
}

impl Simulation {
    /// Returns a structure ready to be simulated, given the `path` to a .json describing both the simulation parameters and the cardiovascular network.
    pub fn new(path: &str) -> Simulation {
        let mut parse = SimulationParsing::read_json(&path.to_string());
        info!("Unmarshalled {}", path);

        let inlet_data = SimulationParsing::read_inlet(&parse.inlet_path);
        info!("Unmarshalled {}", parse.inlet_path);

        let total_time = parse.total_time;

        let consts = Constants::new(&parse);

        info!("Precomputed constants");
        info!(
            ">dt = {:.4e}, c = {:.4e}, cs = {:.4e}, cs2 = {:.4e}",
            consts.dt, consts.c, consts.cs, consts.cs2
        );
        info!(
            ">nu = {:.4e}, tau = {:.4e}, omega = {:.4e}",
            consts.nu, consts.tau, consts.omega
        );

        let mut vessels: Vec<Vessel> = parse
            .vessels
            .iter()
            .map(|v| Vessel::new(&v, consts.clone()))
            .collect();

        for v in &mut vessels {
            for i in 0..v.x_last {
                v.compute_F(i);
            }

            for i in 0..v.x_last {
                v.cells.u[i] += ((consts.dt / 2.0) * v.cells.F[i]) / v.cells.A[i];
                let feq = Simulation::computeFEQ(v.cells.A[i], v.cells.u[i], consts.c);

                v.cells.f0[i] = feq.0;
                v.cells.f1[i] = feq.1;
                v.cells.f2[i] = feq.2;
                v.cells.A[i] = v.cells.f0[i] + v.cells.f1[i] + v.cells.f2[i];
            }
        }

        info!("Initialized vessels");

        Simulation {
            consts,
            total_time,
            vessels,
            inlet: Inlet::new(&inlet_data),
        }
    }

    pub fn computeFEQ(A: f64, u: f64, c: f64) -> (f64, f64, f64) {
        let uc = u / c;
        let uc2 = uc * uc;

        (
            (1.0 / 3.0) * A * (2.0 - 3.0 * uc2),
            (1.0 / 6.0) * A * (1.0 - 3.0 * uc + 3.0 * uc2),
            (1.0 / 6.0) * A * (1.0 + 3.0 * uc + 3.0 * uc2),
        )
    }
}
