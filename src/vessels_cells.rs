//!

use serde::{Deserialize, Serialize};
use std::collections::VecDeque;
use tracing::{event, info, instrument, span, warn, Level};

// #[derive(Debug, Clone)]
// pub struct Populations(pub f64, pub f64, pub f64);

pub type Populations = (f64, f64, f64);
pub type Pop = (VecDeque<f64>, VecDeque<f64>, VecDeque<f64>);

#[derive(Debug, Clone)]
pub struct Cells {
    pub fff: Vec<Populations>,
    /// Stationary populations of particules
    pub f0: VecDeque<f64>,
    /// Right to left populations of particules
    pub f1: VecDeque<f64>,
    /// Left to right populations of particules
    pub f2: VecDeque<f64>,

    /// Area (rho)
    pub A: Vec<f64>,
    /// Velocity
    pub u: Vec<f64>,
    pub stress: Vec<f64>,
    pub deriv: Vec<f64>,

    /// Force
    pub F: Vec<f64>,

    /// Initial area
    pub A0: Vec<f64>,
    /// Squared inversed initial area
    pub s_invA0: Vec<f64>,

    /// Beta
    pub beta: Vec<f64>,

    /// Gamma
    pub gamma: Vec<f64>,
}

impl Default for Cells {
    fn default() -> Self {
        Cells::new(1)
    }
}

impl Cells {
    pub fn new(x_dim: usize) -> Cells {
        Cells {
            fff: vec![(0.0f64, 0.0f64, 0.0f64); x_dim],
            f0: VecDeque::from(vec![0.0f64; x_dim]),
            f1: VecDeque::from(vec![0.0f64; x_dim]),
            f2: VecDeque::from(vec![0.0f64; x_dim]),

            A: vec![0.0f64; x_dim],
            u: vec![0.0f64; x_dim],
            stress: vec![0.0f64; x_dim],
            deriv: vec![0.0f64; x_dim],

            A0: vec![0.0f64; x_dim],
            s_invA0: vec![0.0f64; x_dim],

            F: vec![0.0f64; x_dim],

            beta: vec![0.0f64; x_dim],
            gamma: vec![0.0f64; x_dim],
        }
    }
}
