//! Initializes and runs a simulation

use crate::boundary_conditions::Inlet;
use crate::constants::Constants;
use crate::lbm_algorithm::*;
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
pub struct Simulation<T: Algorithm> {
    /// [Constants](crate::Constants)
    pub consts: Constants,
    /// Duration of the simulation
    pub total_time: f64,
    /// [Vessels](crate::Vessel)
    pub vessels: Vec<Vessel>,
    /// Uniformly sampled inflow (w.r.t. dt)
    pub inlet: Inlet,
    pub algo: T,
}

impl<T: Algorithm> Simulation<T> {
    /// Returns a structure ready to be simulated, given the `path` to a .json describing both the simulation parameters and the cardiovascular network.
    pub fn new(path: &str) -> Simulation<T> {
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

        let mut vessels: Vec<Vessel> = parse.vessels.iter().map(|v| Vessel::new(&v, consts.clone())).collect();

        let algo = T::new(consts.clone());

        for v in &mut vessels {
            v.cells.F = algo.compute_forcing_term(v);
            v.cells.u = algo.compute_velocity(v);
            v.cells.f = algo.compute_FEQ(v);
            v.cells.A = algo.compute_area(v);
        }

        info!("Initialized vessels");

        Simulation {
            consts,
            total_time,
            vessels,
            algo,
            inlet: Inlet::new(&inlet_data),
        }
    }
}
