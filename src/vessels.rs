//! Initializes a vessels

use crate::vessels_cells::Cells;
use serde::{Deserialize, Serialize};
use std::f64::consts::{FRAC_1_PI, PI};
use tracing::{event, info, instrument, span, warn, Level};

/// Represents a vessel
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Vessel {
    /// `id`-th vessel in the vascular network
    pub id: i64,
    /// Name vessel in the vascular network
    pub name: String,
    /// Inflow boundary condition
    #[serde(rename = "root")]
    pub is_inlet: bool,
    /// Length of the vessel [m]
    pub length: f64,
    /// Inlet radius [m]
    pub radius_proximal: f64,
    /// Outlet radius [m]
    pub radius_distal: f64,
    /// Wall thickness of the vessel [m]
    pub wall_thickness: f64,
    /// Young's modulus [Pa]
    pub young_modulus: f64,
    /// Vector of childrens' ids
    pub children: Vec<i64>,
    /// Kind of outflow (`None` if the vessel has children)
    pub outflow: Option<Outflow>,
    #[serde(skip)]
    pub dx: f64,
    #[serde(skip)]
    pub dt: f64,
    #[serde(skip)]
    pub rho: f64,
    #[serde(skip)]
    pub cells: Cells,
    #[serde(skip)]
    pub x_dim: usize,
    #[serde(skip)]
    pub x_last: usize,
}

/// Kind of outflow
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Outflow {
    /// Non-reflective i.e. "as if the vessel goes on infinitely"
    #[serde(rename = "non_reflective")]
    NonReflective,
    /// WK3-lumped parameter model
    #[serde(rename = "wk3")]
    WK3(Wk3),
}

/// Three-element Windkessel model
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Wk3 {
    /// Resistance of blood due to blood viscosity
    pub r1: f64,
    /// Blood inertance
    pub r2: f64,
    /// Compliance of the artery
    pub c: f64,
}

impl Vessel {
    pub fn init(&mut self, dx: f64, dt: f64, rho: f64) {
        self.dx = dx;
        self.dt = dt;
        self.rho = rho;

        self.x_dim = (self.length / dx).round() as usize;
        self.x_last = self.x_dim - 1;

        self.cells = Cells::new(self.x_dim);

        let x_last_f = self.x_last as f64;

        for i in 0..self.x_last {
            let i_f = i as f64;

            // Linear interpolation
            let radius = ((x_last_f - i_f) / (x_last_f)) * self.radius_proximal
                + (i_f / x_last_f) * self.radius_distal;

            self.cells.A[i] = PI * radius * radius;
            self.cells.u[i] = 0.0;

            self.cells.f0[i] = (4.0 / 6.0) * self.cells.A[i];
            self.cells.f1[i] = (1.0 / 6.0) * self.cells.A[i];
            self.cells.f2[i] = (1.0 / 6.0) * self.cells.A[i];

            self.cells.A0[i] = self.cells.A[i];
            self.cells.s_invA0[i] = FRAC_1_PI.sqrt() / radius;
            self.cells.beta[i] = (1.0 / (3.0 * rho * PI.sqrt())) / radius;
        }
    }

    pub fn compute_F(&mut self, cs2: f64, i: usize) {
        let A_derivative = 0.0;

        let first = self.cells.A[i] / self.rho;
        let P_derivative =
            (self.cells.beta[i]) / (2.0 * (self.cells.A[i] * self.cells.A0[i]).sqrt());

        let H_derivative = (self.cells.A[i] / self.rho) * P_derivative;
        let A_derivative = match i {
            // Derivative with i and i+1
            0 => (self.cells.A[i + 1] - self.cells.A[i]) / self.dx,

            // Derivative with i-1 and i-1
            x_last => (self.cells.A[i] - self.cells.A[i - 1]) / self.dx,

            // Derivative with i-1, i and i+1
            _ => (self.cells.A[i + 1] - self.cells.A[i - 1]) / (2.0 * self.dx),
        };

        self.cells.F[i] = A_derivative * ((cs2) - (H_derivative));
    }

    pub fn compute_all_F(&mut self, cs2: f64) {
        (0..self.x_last).for_each(|i| self.compute_F(cs2, i));
    }
}
