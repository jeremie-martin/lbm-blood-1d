use zip_map;
use zip_map_enumerate;

use crate::constants::Constants;
use crate::vessels::Vessel;
use crate::vessels_parsing::*;
use std::collections::VecDeque;
use tracing::{event, info, instrument, span, warn, Level};

/// Defines the LBM steps along a vessel
pub trait Compute {
    /// Constructor
    fn new(consts: Constants) -> Self;

    /// Computes the equilibrium populations
    fn FEQ(&self, A: &f64, u: &f64) -> (f64, f64, f64);
    // fn FEQ(&self, v: &Vessel) -> VecDeque<(f64, f64, f64)>;

    fn BGK(&self, A: &f64, u: &f64, F: &f64, f: (&f64, &f64, &f64)) -> (f64, f64, f64);

    /// Computes the initial macroscopic velocity, without the populations
    fn populations_init(&self, A: &f64, u: &f64, F: &f64, ud: &f64) -> (f64, f64, f64);

    /// Computes the initial macroscopic velocity, without the populations
    fn velocity_init(&self, A: &f64, u: &f64, F: &f64) -> f64;

    /// Computes the macroscopic velocity
    // fn velocity(&self, v: &mut Vessel);
    fn velocity(&self, A: &f64, F: &f64, f: (&f64, &f64, &f64)) -> f64;

    fn stress(&self, A: &f64, u: &f64, F: &f64, beta: &f64, A0: &f64, f: (&f64, &f64, &f64)) -> f64;
    fn dev(&self, v: &mut Vessel);

    /// Computes the macroscopic area
    fn area(&self, f: (&f64, &f64, &f64)) -> f64;

    /// Computes the forcing term given the current macroscopic quantities
    fn forcing_term(&self, v: &mut Vessel, A_lhs: Option<f64>, A_rhs: Option<f64>);

    fn non_reflective<'a>(
        &self,
        v: &'a Vessel,
        A_apply: f64,
        u_apply: f64,
        idx_apply: usize,
        alpha: f64,
    ) -> (f64, f64, f64);

    // /// WK3
    fn wk3(&self, v: &Vessel, wk3: &Wk3) -> (f64, f64);
}

pub struct AlgoBase {
    consts: Constants,
}

impl Compute for AlgoBase {
    fn new(consts: Constants) -> AlgoBase {
        AlgoBase { consts }
    }

    fn FEQ(&self, A: &f64, u: &f64) -> (f64, f64, f64) {
        let uc = u / self.consts.c;
        let uc2 = uc * uc;
        (
            (1.0 / 3.0) * A * (2.0 - 3.0 * uc2),
            (1.0 / 6.0) * A * (1.0 - 3.0 * uc + 3.0 * uc2),
            (1.0 / 6.0) * A * (1.0 + 3.0 * uc + 3.0 * uc2),
        )
    }

    // fn FEQ(&self, v: &Vessel) -> VecDeque<(f64, f64, f64)> {
    //     zip_for_each!(v.cells, |(A, u, mut f)| {
    //         let uc = u / ;
    //         let uc2 = uc * uc;
    //
    //         f.0 = (1.0 / 3.0) * A * (2.0 - 3.0 * uc2);
    //         f.1 = (1.0 / 6.0) * A * (1.0 - 3.0 * uc + 3.0 * uc2);
    //         f.2 = (1.0 / 6.0) * A * (1.0 + 3.0 * uc + 3.0 * uc2);
    //     });
    // }

    fn BGK(&self, A: &f64, u: &f64, F: &f64, f: (&f64, &f64, &f64)) -> (f64, f64, f64) {
        let feq = self.FEQ(A, u);
        // let up = self.consts.c * (f.2 - f.1) / A;
        // let feqp = self.FEQ(A, &up);
        // let feqpd = self.FEQ(A, &(up + (F * self.consts.dt / A)));
        let c = self.consts.c;

        // (
        //     f.0 + (self.consts.omega * (feq.0 - f.0)) + (feqpd.0 - feqp.0) * self.consts.dt,
        //     f.1 + (self.consts.omega * (feq.1 - f.1)) + (feqpd.1 - feqp.1) * self.consts.dt,
        //     f.2 + (self.consts.omega * (feq.2 - f.2)) + (feqpd.2 - feqp.2) * self.consts.dt,
        // )
        (
            f.0 + (self.consts.omega * (feq.0 - f.0)) + (self.consts.omega_force * ((-2.0 * u) / (c * c)) * F),
            f.1 + (self.consts.omega * (feq.1 - f.1)) + (self.consts.omega_force * ((-0.5 * c + u) / (c * c)) * F),
            f.2 + (self.consts.omega * (feq.2 - f.2)) + (self.consts.omega_force * ((0.5 * c + u) / (c * c)) * F),
        )
    }

    fn populations_init(&self, A: &f64, u: &f64, F: &f64, ud: &f64) -> (f64, f64, f64) {
        let feq = self.FEQ(A, u);
        let fneq0 = (4.0 / 6.0) * self.consts.tau_bar * A * ud
            - ((4.0 / 6.0) * self.consts.dt / (2.0 * self.consts.cs2)) * (-0.5 * (2.0 * u * F));

        let fneq1 = (1.0 / 6.0) * (self.consts.tau_bar / self.consts.cs2) * A * (2.0 * self.consts.cs2) * ud
            - ((1.0 / 6.0) * self.consts.dt / (2.0 * self.consts.cs2)) * (self.consts.c * F + (2.0 * u * F));

        let fneq2 = (1.0 / 6.0) * (self.consts.tau_bar / self.consts.cs2) * A * (2.0 * self.consts.cs2) * ud
            - ((1.0 / 6.0) * self.consts.dt / (2.0 * self.consts.cs2)) * (-self.consts.c * F + (2.0 * u * F));

        (feq.0 + fneq0, feq.1 + fneq1, feq.2 + fneq2)
    }

    // TODO: `u +` necessary?
    fn velocity_init(&self, A: &f64, u: &f64, F: &f64) -> f64 {
        u - (F * self.consts.dt) / (2.0 * A)
    }

    fn velocity(&self, A: &f64, F: &f64, f: (&f64, &f64, &f64)) -> f64 {
        ((self.consts.c * (f.2 - f.1)) + (0.5 * F * self.consts.dt)) / A
    }

    fn stress(&self, A: &f64, u: &f64, F: &f64, beta: &f64, A0: &f64, f: (&f64, &f64, &f64)) -> f64 {
        let X = self.consts.dt / (2.0 * self.consts.tau_bar);
        let Y = (1.0 - X);
        let feq = self.FEQ(A, u);
        let bgk = self.BGK(A, u, F, f);

        let first = Y * self.consts.c2 * (f.1 + f.2);
        let second = X * self.consts.c2 * (feq.1 + feq.2);
        let third = X * Y * 2.0 * F * u;

        // let stress =
        //     -(1.0 - (self.consts.dt / (2.0 * self.consts.tau_bar))) * self.consts.c2 * ((f.1 - feq.1) + (f.2 - feq.2));
        // let stress = self.consts.c2 * ((f.1 - feq.1) + (f.2 - feq.2));
        // let ud = -stress / (A * self.consts.cs2 * self.consts.tau_bar);
        let first = (-Y) * self.consts.c2 * ((f.1 - feq.1) + (f.2 - feq.2));
        let second = (-(0.5 * self.consts.dt) * Y) * 2.0 * F * u;
        let stress = first + second;
        let nu = (self.consts.tau_bar - (self.consts.dt / 2.0)) * self.consts.cs2;
        let ud = stress / (A * nu * 2.0);

        let dP = A * ud;

        // (2.0 * A0 * dP * (A / A0).sqrt()) / (beta)

        let dUA = -((self.consts.c * (f.2 - f.1)) + (0.5 * F * self.consts.dt)) / (A * A);

        let ddd = (self.consts.c * (f.0 + 2.0 * (f.1)) - self.consts.dt * 0.5 * F)
            / ((self.consts.c - u) * (self.consts.c - u));

        ud
        //
    }
    // fn velocity(&self, v: &mut Vessel) {
    //     zip_for_each!(v.cells, |(A, mut u, F, f)| {
    //         *u = self._velocity(A, F, f);
    //         // (self.consts.c * (f.2 - f.1) + ((self.consts.dt / 3.0) * F)) / A
    //     });
    // }

    fn area(&self, f: (&f64, &f64, &f64)) -> f64 {
        return f.0 + f.1 + f.2;
    }
    fn dev(&self, v: &mut Vessel) {
        for i in 0..v.x_dim {
            let A_derivative = match i {
                0 => {
                    v.cells.u[i + 1] - v.cells.u[i]

                    // (-25.0 / 12.0) * v.cells.A[i] + 4.0 * v.cells.A[i + 1] - 3.0 * v.cells.A[i + 2]
                    //     + (4.0 / 3.0) * v.cells.A[i + 3]
                    //     - (1.0 / 4.0) * v.cells.A[i + 4]
                }
                i if (i == v.x_last) => {
                    v.cells.u[i] - v.cells.u[i - 1]
                    // (25.0 / 12.0) * v.cells.A[i] - 4.0 * v.cells.A[i - 1] + 3.0 * v.cells.A[i - 2]
                    //     - (4.0 / 3.0) * v.cells.A[i - 3]
                    //     + (1.0 / 4.0) * v.cells.A[i - 4]
                }
                1 => -0.5 * v.cells.u[i - 1] + 0.5 * v.cells.u[i + 1],
                i if (i == v.x_last - 1) => -0.5 * v.cells.u[i - 1] + 0.5 * v.cells.u[i + 1],
                _ => {
                    // -0.5 * v.cells.A[i - 1] + 0.5 * v.cells.A[i + 1]
                    (1.0 / 12.0) * v.cells.u[i - 2] - (2.0 / 3.0) * v.cells.u[i - 1] + (2.0 / 3.0) * v.cells.u[i + 1]
                        - (1.0 / 12.0) * v.cells.u[i + 2]
                }
            };

            let A_derivative = A_derivative / (self.consts.dx);
            v.cells.deriv[i] = A_derivative;
        }
    }

    // TODO: make a A_derivative function
    fn forcing_term(&self, v: &mut Vessel, A_lhs: Option<f64>, A_rhs: Option<f64>) {
        for i in 0..v.x_dim {
            let A_derivative = match i {
                0 if A_lhs.is_none() => {
                    v.cells.A[i + 1] - v.cells.A[i]

                    // (-25.0 / 12.0) * v.cells.A[i] + 4.0 * v.cells.A[i + 1] - 3.0 * v.cells.A[i + 2]
                    //     + (4.0 / 3.0) * v.cells.A[i + 3]
                    //     - (1.0 / 4.0) * v.cells.A[i + 4]
                }
                i if (i == v.x_last && A_rhs.is_none()) => {
                    v.cells.A[i] - v.cells.A[i - 1]
                    // 1.5 * v.cells.A[i] - 2.0 * v.cells.A[i - 1] + 0.5 * v.cells.A[i - 2]
                    // (25.0 / 12.0) * v.cells.A[i] - 4.0 * v.cells.A[i - 1] + 3.0 * v.cells.A[i - 2]
                    //     - (4.0 / 3.0) * v.cells.A[i - 3]
                    //     + (1.0 / 4.0) * v.cells.A[i - 4]
                }
                0 if A_lhs.is_some() => -0.5 * A_lhs.unwrap() + 0.5 * v.cells.A[i + 1],
                i if (i == v.x_last && A_rhs.is_some()) => -0.5 * v.cells.A[i - 1] + 0.5 * A_rhs.unwrap(),
                _ => -0.5 * v.cells.A[i - 1] + 0.5 * v.cells.A[i + 1],
                // 1 => -0.5 * v.cells.A[i - 1] + 0.5 * v.cells.A[i + 1],
                // i if (i == v.x_last - 1) => -0.5 * v.cells.A[i - 1] + 0.5 * v.cells.A[i + 1],
                // _ => {
                //     // -0.5 * v.cells.A[i - 1] + 0.5 * v.cells.A[i + 1]
                //     (1.0 / 12.0) * v.cells.A[i - 2] - (2.0 / 3.0) * v.cells.A[i - 1] + (2.0 / 3.0) * v.cells.A[i + 1]
                //         - (1.0 / 12.0) * v.cells.A[i + 2]
                // }
            };

            let A_derivative = A_derivative / (self.consts.dx);

            v.cells.F[i] = -A_derivative
                * (((v.cells.beta[i] * (v.cells.A[i] / v.cells.A0[i]).sqrt()) / (2.0 * self.consts.rho))
                    - self.consts.cs2);
            // v.cells.F[i] = 0.5
            //     * (2.0 * self.consts.cs2 * self.consts.rho - v.cells.beta[i] * (v.cells.A[i] / v.cells.A0[i]).sqrt())
            //     * A_derivative
            //     / self.consts.rho;
        }

        // zip_map_enumerate!(v.cells, |(i, A0, beta)| {
        //     let P_derivative = (beta) / (2.0 * (v.cells.A[i] * A0).sqrt());
        //     let H_derivative = (v.cells.A[i] / self.consts.rho) * P_derivative;
        //
        //     let A_derivative = match i {
        //         0 => (v.cells.A[i + 1] - v.cells.A[i]) / self.consts.dx,
        //         i if (i == v.x_last) => (v.cells.A[i] - v.cells.A[i - 1]) / self.consts.dx,
        //         _ => (v.cells.A[i + 1] - v.cells.A[i - 1]) / (2.0 * self.consts.dx),
        //     };
        //
        //     A_derivative * ((self.consts.cs2) - (H_derivative))
        // })
    }

    fn non_reflective<'a>(
        &self,
        v: &'a Vessel,
        A_apply: f64,
        u_apply: f64,
        idx_apply: usize,
        alpha: f64,
    ) -> (f64, f64, f64) {
        let u = u_apply;
        // let u = self.velocity(
        //     &v.cells.A[v.x_last - 1],
        //     &v.cells.F[v.x_last - 1],
        //     (
        //         &v.cells.f0[v.x_last - 1],
        //         &v.cells.f1[v.x_last - 1],
        //         &v.cells.f2[v.x_last - 1],
        //     ),
        // );

        let feq = self.FEQ(&v.cells.A[idx_apply], &v.cells.u[idx_apply]);
        let feqn = self.FEQ(&A_apply, &u);

        // let f1 = if (idx_apply == v.x_last && false) {
        //     ((self.consts.dt * 0.5 / self.consts.c) * v.cells.F[idx_apply]) + v.cells.f2[idx_apply]
        //         - A_apply * (u / self.consts.c)
        // } else {
        //     v.cells.f1[idx_apply]
        // };

        let f0 = v.cells.f0[idx_apply];
        let f1 = v.cells.f1[idx_apply];
        let f2 = v.cells.f2[idx_apply];

        (
            (f0 - feq.0 + feqn.0) * alpha + f0 * (1.0 - alpha),
            (f1 - feq.1 + feqn.1) * alpha + f1 * (1.0 - alpha),
            (f2 - feq.2 + feqn.2) * alpha + f2 * (1.0 - alpha),
        )
    }

    fn wk3(&self, v: &Vessel, wk3: &Wk3) -> (f64, f64) {
        (0.0, 0.0)
        // let P = v.cells.beta[v.x_last] * ((v.cells.A[v.x_last] / v.cells.A0[v.x_last]).sqrt() - 1.0);
        // let Q = v.cells.A[v.x_last] * v.cells.u[v.x_last];
        //
        // let Qdt = ((3.0 / 2.0) + v.Q_old[0] - 2.0 * v.Q_old[1] + 0.5 * v.Q_old[2]) / self.consts.dt;
        //
        // let P_new =
        //     P + self.consts.dt * (((1.0 + wk3.r1 / wk3.r2) * Q) + (wk3.c * wk3.r1 * Qdt) - (P / wk3.r2)) / wk3.c;
        //
        // (P_new, Q)

        // let aL = 2.0 * v.cells.F[v.x_last - 1] - v.cells.F[v.x_last - 2];
        // let f0 = v.cells.f0[v.x_last];
        // let f1 = v.cells.f1[v.x_last];
        // let f2 = v.cells.f2[v.x_last];
        // let A0 = v.cells.A0[v.x_last];
        // let A = v.cells.A[v.x_last];
        // let beta = v.cells.beta[v.x_last];
        // let rho = v.cells.rho[v.x_last];
        // let cpulse = 3.0 * v.cells.gamma[v.x_last] * A.sqrt() * 0.5;
        // let df1 = f1 - (1.0 / 6.0) * A0;
        // let df2 = f2 - (1.0 / 6.0) * A0;
        //
        // let idx = v.x_last - ((cpulse * self.consts.dt) as usize);
        // let pressure = v.cells.f2[idx] * ((v.cells.A[idx] / v.cells.A0[idx]).sqrt() - 1.0);
        // let flow = v.cells.A[idx] * v.cells.u[idx];
        // let pf = 0.5 * (pressure + wk3.r1 * flow);
        // let pb = wk3.r2 * self.consts.dt * pf + wk3.r1 * wk3.r2 * wk3.c * (pressure - wk3.r1 * flow);
        // let pb = pb / ((wk3.r1 + wk3.r2) * self.consts.dt + 2.0 * wk3.r1 * wk3.r2 * wk3.c);
        // let C0 = (A0 / (self.consts.rho * cpulse));
        // let q = C0 *
        //
        // let _f0 = (4.0 / 6.0) * A0 - 2.0 * df2 + C0 * (pf + pb)
        //     - (self.consts.dt * aL / 2.0)
        //     + q;
        // let _f1 = f2 + (self.consts.dt * aL / 2.0) - q;

        // (_f1, _f1)
    }
}
