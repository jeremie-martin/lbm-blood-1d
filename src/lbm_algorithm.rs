use zip_map;
use zip_map_enumerate;

use crate::constants::Constants;
use crate::vessels::Vessel;
use crate::vessels_cells::Populations;

/// Defines the LBM steps along a vessel
pub trait Algorithm {
    /// Constructor
    fn new(consts: Constants) -> Self;

    /// Computes the equilibrium populations
    fn compute_FEQ(&self, v: &Vessel) -> Vec<Populations>;

    /// Computes the macroscopic velocity
    fn compute_velocity(&self, v: &Vessel) -> Vec<f64>;

    /// Computes the macroscopic area
    fn compute_area(&self, v: &Vessel) -> Vec<f64>;

    /// Computes the forcing term given the current macroscopic quantities
    fn compute_forcing_term(&self, v: &Vessel) -> Vec<f64>;
}

pub struct AlgoBase {
    consts: Constants,
}

impl Algorithm for AlgoBase {
    fn new(consts: Constants) -> AlgoBase {
        AlgoBase { consts }
    }

    fn compute_FEQ(&self, v: &Vessel) -> Vec<Populations> {
        zip_map!(v.cells, |(A, u)| {
            let uc = u / self.consts.c;
            let uc2 = uc * uc;
            (
                (1.0 / 3.0) * A * (2.0 - 3.0 * uc2),
                (1.0 / 6.0) * A * (1.0 - 3.0 * uc + 3.0 * uc2),
                (1.0 / 6.0) * A * (1.0 + 3.0 * uc + 3.0 * uc2),
            )
        })
    }

    fn compute_velocity(&self, v: &Vessel) -> Vec<f64> {
        zip_map!(v.cells, |(A, u, F)| u + ((self.consts.dt / 2.0) * F) / A)
    }

    fn compute_area(&self, v: &Vessel) -> Vec<f64> {
        zip_map!(v.cells, |f| f.0 + f.1 + f.2)
    }

    fn compute_forcing_term(&self, v: &Vessel) -> Vec<f64> {
        zip_map_enumerate!(v.cells, |(i, A0, beta)| {
            let P_derivative = (beta) / (2.0 * (v.cells.A[i] * A0).sqrt());
            let H_derivative = (v.cells.A[i] / self.consts.rho) * P_derivative;

            let A_derivative = match i {
                0 => (v.cells.A[i + 1] - v.cells.A[i]) / self.consts.dx,
                i if (i == v.x_last) => (v.cells.A[i] - v.cells.A[i - 1]) / self.consts.dx,
                _ => (v.cells.A[i + 1] - v.cells.A[i - 1]) / (2.0 * self.consts.dx),
            };

            A_derivative * ((self.consts.cs2) - (H_derivative))
        })
    }
}
