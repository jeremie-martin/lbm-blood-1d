use serde::{Deserialize, Serialize};
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
}
