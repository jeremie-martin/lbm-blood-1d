//! Initializes a vessels

use crate::constants::Constants;
use crate::vessels_cells::*;
use crate::vessels_parsing::*;
use serde::{Deserialize, Serialize};
use slotmap::{new_key_type, Key, SlotMap, Slottable};
use std::f64::consts::{FRAC_1_PI, PI};
use tracing::{event, info, instrument, span, warn, Level};

new_key_type! {
    pub struct VesselKey;
}

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
    pub children: Vec<usize>,
    /// Kind of outflow (`None` if the vessel has children)
    pub outflow: Option<Outflow>,
    pub consts: Constants,
    pub lhs_give: VesselKey,
    pub rhs_give: VesselKey,
    pub lhs_recv: Vec<VesselKey>,
    pub rhs_recv: Vec<VesselKey>,
    pub cells: Cells,
    pub x_dim: usize,
    pub x_last: usize,
}

impl Vessel {
    pub fn new(parse: &VesselParsing, consts: Constants) -> Vessel {
        let x_dim = (parse.length / consts.dx).round() as usize;
        let x_last = x_dim - 1;

        let mut outlet = parse.outflow.clone();

        match outlet {
            Some(Outflow::WK3(ref mut wk3)) => {
                wk3.pc = 0.0;
            }
            _ => (),
        }

        let mut cells = Cells::new(x_dim);

        let x_last_f = x_last as f64;
        let slope = (parse.radius_distal - parse.radius_proximal) / parse.length;

        let ah = 0.2802;
        let bh = -5.053e2;
        let ch = 0.1324;
        let dh = -0.1114e2;

        let Rm = (parse.radius_distal + parse.radius_proximal) * 0.5;

        let h0 = parse.wall_thickness;
        info!("h0 before {}", h0);
        let h0 = Rm * (ah * (bh * Rm).exp() + ch * (dh * Rm).exp());
        info!("h0 after {}", h0);

        for i in 0..x_dim {
            let i_f = i as f64;

            // Linear interpolation
            let radius = slope * i_f * consts.dx + parse.radius_proximal;

            cells.A[i] = PI * radius * radius;
            cells.u[i] = 0.0;

            cells.f0[i] = (4.0 / 6.0) * cells.A[i];
            cells.f1[i] = (1.0 / 6.0) * cells.A[i];
            cells.f2[i] = (1.0 / 6.0) * cells.A[i];

            cells.A0[i] = cells.A[i];
            cells.s_invA0[i] = (1.0 / cells.A0[i]).sqrt();
            cells.beta[i] = cells.s_invA0[i] * h0 * PI.sqrt() * parse.young_modulus / 0.75;
            cells.gamma[i] = cells.beta[i] * (1.0 / (3.0 * consts.rho * PI.sqrt())) / radius;
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
            children: parse.children.clone().iter().map(|id| id - 1).collect(),
            lhs_give: VesselKey::null(),
            rhs_give: VesselKey::null(),
            lhs_recv: Vec::<VesselKey>::new(),
            rhs_recv: Vec::<VesselKey>::new(),
            consts,
            outflow: outlet,
            x_dim,
            x_last,
            cells,
        }
    }
}
