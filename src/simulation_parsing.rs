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
pub struct SimulationParsing {
    /// Space step [m]
    pub dx: f64,
    /// Wall visco-elasticity (Pa.s.m-1)
    pub gamma: f64,
    /// Blood dynamic viscosity [Pa.s = kg.m-1.s-1]
    pub mu: f64,
    /// Blood density [kg.m-3]
    pub rho: f64,
    /// Time of the simulation
    pub total_time: f64,
    /// Path to the inlet file
    pub inlet_path: String,
    /// [Vessels](crate::Vessel)
    pub vessels: Vec<Vessel>,
}

impl SimulationParsing {
    /// Unmarshalls a json file describing both the simulation parameters and the cardiovascular network.
    pub fn read_json(path: &String) -> SimulationParsing {
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
}
