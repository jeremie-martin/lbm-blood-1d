//! Initializes a vessels

use crate::constants::Constants;
use crate::vessels_cells::*;
use crate::vessels_parsing::*;
use serde::{Deserialize, Serialize};
use std::f64::consts::{FRAC_1_PI, PI};
use tracing::{event, info, instrument, span, warn, Level};

/// Represents a vessel
#[derive(Debug, Clone)]
pub struct Vessel {
    /// `id`-th vessel in the vascular network
    pub id: i64,
    /// Name vessel in the vascular network
    pub name: String,
    /// Inflow boundary condition
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
    pub consts: Constants,
    pub cells: Cells,
    pub x_dim: usize,
    pub x_last: usize,
}

impl Vessel {
    pub fn new(parse: &VesselParsing, consts: Constants) -> Vessel {
        let x_dim = (parse.length / consts.dx).round() as usize;
        let x_last = x_dim - 1;

        let mut cells = Cells::new(x_dim);

        let x_last_f = x_last as f64;

        for i in 0..x_last {
            let i_f = i as f64;

            // Linear interpolation
            let radius =
                ((x_last_f - i_f) / (x_last_f)) * parse.radius_proximal + (i_f / x_last_f) * parse.radius_distal;

            cells.A[i] = PI * radius * radius;
            cells.u[i] = 0.0;

            cells.f[i].0 = (4.0 / 6.0) * cells.A[i];
            cells.f[i].1 = (1.0 / 6.0) * cells.A[i];
            cells.f[i].2 = (1.0 / 6.0) * cells.A[i];

            cells.f0[i] = (4.0 / 6.0) * cells.A[i];
            cells.f1[i] = (1.0 / 6.0) * cells.A[i];
            cells.f2[i] = (1.0 / 6.0) * cells.A[i];

            cells.A0[i] = cells.A[i];
            cells.s_invA0[i] = FRAC_1_PI.sqrt() / radius;
            cells.beta[i] = (1.0 / (3.0 * consts.rho * PI.sqrt())) / radius;
        }

        Vessel {
            id: parse.id - 1,
            name: parse.name.clone(),
            is_inlet: parse.is_inlet,
            length: parse.length,
            radius_proximal: parse.radius_proximal,
            radius_distal: parse.radius_distal,
            wall_thickness: parse.wall_thickness,
            young_modulus: parse.young_modulus,
            children: parse.children.clone(),
            outflow: parse.outflow.clone(),
            consts,
            x_dim,
            x_last,
            cells,
        }
    }

    pub fn computeFEQ(&self, A: f64, u: f64) -> Populations {
        let uc = u / self.consts.c;
        let uc2 = uc * uc;

        (
            (1.0 / 3.0) * A * (2.0 - 3.0 * uc2),
            (1.0 / 6.0) * A * (1.0 - 3.0 * uc + 3.0 * uc2),
            (1.0 / 6.0) * A * (1.0 + 3.0 * uc + 3.0 * uc2),
        )
    }

    pub fn computeEQ(&self, i: usize) -> Populations {
        let uc = self.cells.u[i] / self.consts.c;
        let uc2 = uc * uc;

        (
            (1.0 / 3.0) * self.cells.A[i] * (2.0 - 3.0 * uc2),
            (1.0 / 6.0) * self.cells.A[i] * (1.0 - 3.0 * uc + 3.0 * uc2),
            (1.0 / 6.0) * self.cells.A[i] * (1.0 + 3.0 * uc + 3.0 * uc2),
        )
    }

    pub fn update_velocity(&self, i: usize) -> f64 {
        self.cells.u[i] + ((self.consts.dt / 2.0) * self.cells.F[i]) / self.cells.A[i]
    }

    pub fn update_area(&self, i: usize) -> f64 {
        self.cells.f[i].0 + self.cells.f[i].1 + self.cells.f[i].2
    }

    pub fn populations_init(&mut self, i: usize) {
        self.cells.u[i] += ((self.consts.dt / 2.0) * self.cells.F[i]) / self.cells.A[i];

        self.cells.f[i] = self.computeFEQ(self.cells.A[i], self.cells.u[i]);

        // self.cells.f0[i] = feq.0;
        // self.cells.f1[i] = feq.1;
        // self.cells.f2[i] = feq.2;
        self.cells.A[i] = self.cells.f[i].0 + self.cells.f[i].1 + self.cells.f[i].2;
    }

    pub fn compute_F2(&self, i: usize) -> f64 {
        let first = self.cells.A[i] / self.consts.rho;
        let P_derivative = (self.cells.beta[i]) / (2.0 * (self.cells.A[i] * self.cells.A0[i]).sqrt());

        let H_derivative = (self.cells.A[i] / self.consts.rho) * P_derivative;

        let A_derivative = match i {
            // Derivative with i and i+1
            0 => (self.cells.A[i + 1] - self.cells.A[i]) / self.consts.dx,

            // Derivative with i-1 and i-1
            i if (i == self.x_last) => (self.cells.A[i] - self.cells.A[i - 1]) / self.consts.dx,

            // Derivative with i-1, i and i+1
            _ => (self.cells.A[i + 1] - self.cells.A[i - 1]) / (2.0 * self.consts.dx),
        };

        A_derivative * ((self.consts.cs2) - (H_derivative))
    }

    pub fn compute_F(&mut self, i: usize) {
        let first = self.cells.A[i] / self.consts.rho;
        let P_derivative = (self.cells.beta[i]) / (2.0 * (self.cells.A[i] * self.cells.A0[i]).sqrt());

        let H_derivative = (self.cells.A[i] / self.consts.rho) * P_derivative;

        let A_derivative = match i {
            // Derivative with i and i+1
            0 => (self.cells.A[i + 1] - self.cells.A[i]) / self.consts.dx,

            // Derivative with i-1 and i-1
            i if (i == self.x_last) => (self.cells.A[i] - self.cells.A[i - 1]) / self.consts.dx,

            // Derivative with i-1, i and i+1
            _ => (self.cells.A[i + 1] - self.cells.A[i - 1]) / (2.0 * self.consts.dx),
        };

        self.cells.F[i] = A_derivative * ((self.consts.cs2) - (H_derivative));
    }

    pub fn compute_all_F(&mut self, cs2: f64) {
        (0..self.x_last).for_each(|i| self.compute_F(i));
    }
}
