//! Initializes and runs a simulation

use crate::vessels::*;
use itertools_num::linspace;
use serde::{Deserialize, Serialize};
use splines::{Interpolation, Key, Spline};
use std::fs::File;
use tracing::{event, info, instrument, span, warn, Level};

/// Represent the couple `(time, volumetric_flow_rate)`
pub type InletRaw = (f64, f64);

/// Contains the simulation parameters and the vessels.
#[derive(Debug, Serialize, Deserialize)]
pub struct Simulation {
    /// Space step [m]
    pub dx: f64,
    /// Time step [m]
    #[serde(skip)]
    pub dt: f64,
    /// Wall visco-elasticity (Pa.s.m^-1)
    pub gamma: f64,
    /// Blood dynamic viscosity [Pa.s = kg.m-1.s-1]
    pub mu: f64,
    /// Blood density [kg.m^-3]
    pub rho: f64,
    /// Time of the simulation
    pub total_time: f64,
    /// [Vessels](crate::Vessel)
    pub vessels: Vec<Vessel>,
    /// Path to the inlet file
    pub inflow_path: String,
    /// Uniformly sampled inflow (w.r.t. dt)
    #[serde(skip)]
    pub inflow: Vec<f64>,
}

impl Simulation {
    /// Returns a structure ready to be simulated, given the `path` to a .json describing both the simulation parameters and the cardiovascular network.
    pub fn new(path: &str) -> Simulation {
        let mut simulation = Simulation::read_json(&path.to_string());
        info!("Unmarshalled {}", path);

        let inlet = Simulation::read_inlet(&simulation.inflow_path);
        info!("Unmarshalled {}", simulation.inflow_path);

        simulation.inflow = simulation.standardize_inlet(&inlet);
        info!("Uniformly interpolated inflow");

        simulation
    }

    /// Unmarshalls a json file describing both the simulation parameters and the cardiovascular network.
    pub fn read_json(path: &String) -> Simulation {
        let json = File::open(path).expect("File not found");
        serde_json::from_reader(&json).unwrap()
    }

    /// Unmarshalls an inlet file
    pub fn read_inlet(path: &String) -> Vec<InletRaw> {
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b' ')
            .from_path(path)
            .unwrap();

        rdr.deserialize().map(|r| r.unwrap()).collect()
    }

    /// Transforms a vector of points `(x, y)` to a vector of evenly spaced `y`, w.r.t. `dt`.
    ///
    /// An inlet file defines the volumetric flow rate of the input heart wave for a given time `T`.
    ///
    /// This method
    ///     - Interpolates the values at `(T / dt).round() - 1` points (the last one is the same as
    ///     the first one at t=0)
    ///     - Updates the time step accordingly: `dt = T / (T / dt).round()`
    ///
    /// Returns a vector of uniformly sampled volumetric flow rate, getting ridded of the axis
    pub fn standardize_inlet(&mut self, inlet: &Vec<InletRaw>) -> Vec<f64> {
        let keys = inlet[1..]
            .iter()
            .map(|(x, y)| Key::new(*x, *y, Interpolation::Linear))
            .collect();

        let spline = Spline::from_vec(keys);
        let T = inlet.last().unwrap().0;
        let N = (T / self.dt).round();
        self.dt = T / N;

        linspace(0.0, T - self.dt, (N - 1.0) as usize)
            .map(|x| spline.clamped_sample(x).unwrap())
            .collect()
    }
}
