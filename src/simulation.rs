//! Initializes and runs a simulation

use crate::boundary_conditions::Inlet;
use crate::simulation_parsing::*;
use crate::vessels::*;
use itertools_num::linspace;
use serde::{Deserialize, Serialize};
use splines::{Interpolation, Key, Spline};
use std::fs::File;
use tracing::{event, info, instrument, span, warn, Level};

/// Contains the simulation parameters and the vessels.
#[derive(Debug, Serialize, Deserialize)]
pub struct Simulation {
    /// Space step [m]
    pub dx: f64,
    /// Time step, depends on dx [m]
    pub dt: f64,
    /// Lattice velocity [m.s-1]
    pub c: f64,
    /// Lattice velocity squared [m2.s-2]
    pub c2: f64,
    /// Lattice speed of sound [m.s-1]
    pub cs: f64,
    /// Lattice speed of sound squared [m2.s-1]
    pub cs2: f64,
    /// Relaxation rate [dimensionless]
    pub omega: f64,
    /// Time of the simulation
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

        let dx = parse.dx;
        let dt = dx * dx;

        let total_time = parse.total_time;

        let mu = parse.mu;
        let rho = parse.rho;

        let c = dx / dt;
        let c2 = c * c;
        let cs = c / 3.0f64.sqrt();
        let cs2 = c2 / 3.0;

        // Kinematic viscosity (ratio of a fluid's dynamic viscosity to the fluid's density) [m2.s-1]
        let nu = mu / rho;

        // Relaxation time [s]
        let tau = nu / cs2;

        let omega = dt / (tau + (dt / 2.0));

        info!("Precomputed constants");
        info!(
            "dt = {:.4e}, c = {:.4e}, cs = {:.4e}, cs2 = {:.4e}, nu = {:.4e}, tau = {:.4e}, omega = {:.4e}",
            dt, c, cs, cs2, nu, tau, omega
        );

        Simulation {
            dx,
            dt,
            c,
            c2,
            cs,
            cs2,
            omega,
            total_time,
            inlet: Inlet::new(&inlet_data),
            vessels: parse.vessels,
        }
    }
}
