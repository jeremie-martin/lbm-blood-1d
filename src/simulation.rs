use crate::vessels::*;
use itertools_num::linspace;
use serde::{Deserialize, Serialize};
use splines::{Interpolation, Key, Spline};
use std::fs::File;
use tracing::{event, info, instrument, span, warn, Level};

#[derive(Debug, Serialize, Deserialize)]
pub struct Simulation {
    pub dx: f64,
    pub gamma: f64,
    pub inflow_path: String,
    #[serde(skip)]
    pub inflow: Vec<f64>,
    pub my: f64,
    pub rho: f64,
    pub total_time: f64,
    pub vessels: Vec<Vessel>,
}

impl Simulation {
    pub fn new(path: &str) -> Simulation {
        let mut simulation = Simulation::read_json(path);
        info!("Unmarshalled {}", path);

        let inlet = Simulation::read_inlet(&simulation.inflow_path);
        info!("Unmarshalled {}", simulation.inflow_path);

        simulation.inflow = Simulation::standardize_inlet(&inlet);
        info!("Uniformly interpolated inflow");

        simulation
    }

    pub fn read_json(path: &str) -> Simulation {
        let json = File::open(path).expect("File not found");
        serde_json::from_reader(&json).unwrap()
    }

    pub fn read_inlet(path: &String) -> Vec<(f64, f64)> {
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b' ')
            .from_path(path)
            .unwrap();

        rdr.deserialize().map(|r| r.unwrap()).collect()
    }

    pub fn standardize_inlet(inlet: &Vec<(f64, f64)>) -> Vec<f64> {
        let keys = inlet
            .iter()
            .map(|(x, y)| Key::new(*x, *y, Interpolation::Linear))
            .collect();

        let spline = Spline::from_vec(keys);
        let low = inlet.first().unwrap().0;
        let high = inlet.last().unwrap().0;

        linspace(low, high, 20)
            .map(|x| spline.clamped_sample(x).unwrap())
            .collect()
    }
}
