use serde::{Deserialize, Serialize};
use std::collections::VecDeque;
use std::fs::File;
use tracing::{event, info, instrument, span, warn, Level};

/// Contains the simulation parameters and the vessels.
#[derive(Debug, Serialize, Deserialize)]
pub struct VesselParsing {
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
    #[serde(default)]
    /// TODO
    pub pc: f64,
    #[serde(default)]
    pub P_old: VecDeque<f64>,
    #[serde(default)]
    pub Q_old: VecDeque<f64>,
}

impl Wk3 {
    pub fn owo(&self, A: f64, u: f64, A0: f64, beta: f64, dt: f64) -> (f64, f64, f64, f64) {
        let P = beta * ((A / A0).sqrt() - 1.0);
        let Q = A * u;

        let Qdt = ((3.0 / 2.0) * Q - 2.0 * self.Q_old[0] + 0.5 * self.Q_old[1]) / dt;
        let Pdt = ((3.0 / 2.0) * P - 2.0 * self.P_old[0] + 0.5 * self.P_old[1]) / dt;

        let P1 = (((1.0 + self.r1 / self.r2) * Q) + (self.r1 * Qdt * self.c) - (P / self.r2)) / self.c;
        let Q_new = Q + dt * ((P / self.r2) + (self.c * Pdt) - ((1.0 + self.r1 / self.r2) * Q)) / (self.c * self.r1);

        let Qdt2 = ((3.0 / 2.0) * Q_new - 2.0 * Q + 0.5 * self.Q_old[0]) / dt;
        let P2 = (((1.0 + self.r1 / self.r2) * Q_new) + (self.r1 * Qdt2 * self.c) - (P / self.r2)) / self.c;
        let P_new = P + 0.5 * dt * (P1 + P2);
        // let P_new = P + dt * (((1.0 + self.r1 / self.r2) * Q) + (self.r1 * Qdt * self.c) - (P / self.r2)) / self.c;

        (P, P_new, Q, Q_new)
    }
    pub fn compute(
        &mut self,
        v_A: f64,
        v_u: f64,
        v_F: f64,
        v_A0: f64,
        v_beta: f64,
        v_gamma: f64,
        dt: f64,
    ) -> (f64, f64) {
        let Al = v_A;
        let ul = v_u;

        let f1 = (1.0 / 6.0) * v_A0;
        self.pc += (dt / self.c) * (Al * ul - self.pc / self.r2);
        let ssAl = Al.sqrt().sqrt();
        let sgamma = 2.0 * (6.0 * v_gamma).sqrt();
        let sA0 = v_A0.sqrt();
        let bA0 = v_beta / sA0;

        let mut xn = Al;

        let mut count = 0;
        loop {
            let x0 = xn;

            let f = x0 * self.r1 * (ul + sgamma * (ssAl - x0.sqrt().sqrt())) - (bA0 * (x0.sqrt() - sA0)) + self.pc;
            let df = self.r1 * (ul + sgamma * (ssAl - 1.25 * x0.sqrt().sqrt())) - bA0 * 0.5 / x0.sqrt();

            xn = x0 - f / df;

            count += 1;
            if (xn - x0).abs() <= 1e-6 {
                break;
            }
        }

        let us = (v_beta * ((xn / v_A0).sqrt() - 1.0)) / (xn * self.r1);
        if count > 1 {
            info!("owo { }{} {} {} {} {}", v_A, xn, v_u, us, count, self.pc);
        }
        (xn, us)
    }
}
