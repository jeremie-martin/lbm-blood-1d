//! This module

use crate::simulation_parsing::*;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Constants {
    /// Space step [m]
    pub dx: f64,
    /// Time step, depends on dx [m]
    pub dt: f64,
    /// Blood dynamic viscosity [Pa.s = kg.m-1.s-1]
    pub mu: f64,
    /// Blood density [kg.m-3]
    pub rho: f64,
    /// Lattice velocity [m.s-1]
    pub c: f64,
    /// Lattice velocity squared [m2.s-2]
    pub c2: f64,
    /// Lattice speed of sound [m.s-1]
    pub cs: f64,
    /// Lattice speed of sound squared [m2.s-1]
    pub cs2: f64,
    /// Kinematic viscosity (ratio of a fluid's dynamic viscosity to the fluid's density) [m2.s-1]
    pub nu: f64,
    /// Relaxation time [s]
    pub tau: f64,
    /// Adjusted Relaxation time [s]
    pub tau_bar: f64,
    /// Relaxation rate [dimensionless]
    pub omega: f64,
    /// Relaxation rate of the Force [dimensionless]
    pub omega_force: f64,
}

impl Constants {
    pub fn new(parse: &SimulationParsing) -> Constants {
        let dx = parse.dx;
        let dt = dx * dx;

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

        // Adjusted tau
        let tau_bar = tau + (dt / 2.0);

        let omega = dt / tau_bar;
        let omega_force = dt * (1.0 - (dt / (2.0 * tau_bar)));

        Constants {
            dx,
            dt,
            mu,
            rho,
            c,
            cs,
            c2,
            cs2,
            nu,
            tau,
            tau_bar,
            omega,
            omega_force,
        }
    }
}
