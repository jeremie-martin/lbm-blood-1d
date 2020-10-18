use serde::{Deserialize, Serialize};
use tracing::{event, info, instrument, span, warn, Level};

#[derive(Debug, Serialize, Deserialize)]
pub struct Simulation {
    pub dx: f64,
    pub gamma: f64,
    pub inflow_path: String,
    pub my: f64,
    pub rho: f64,
    pub total_time: f64,
    pub vessels: Vec<Vessel>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Vessel {
    pub id: i64,
    pub name: String,
    pub root: Option<bool>,
    pub length: f64,
    pub radius_proximal: f64,
    pub radius_distal: f64,
    #[serde(rename = "as")]
    pub _as: f64,
    pub ps: f64,
    pub wall_thickness: f64,
    pub young_modulus: f64,
    pub children: Vec<i64>,
    pub outflow: OutflowUnion,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct OutflowClass {
    pub wk3: Wk3,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Wk3 {
    pub c: f64,
    pub rc: f64,
    pub z: f64,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(untagged)]
pub enum OutflowUnion {
    Enum(OutflowEnum),
    OutflowClass(OutflowClass),
}

#[derive(Debug, Serialize, Deserialize)]
pub enum OutflowEnum {
    #[serde(rename = "non_reflective")]
    NonReflective,
}
