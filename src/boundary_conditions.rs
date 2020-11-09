//! Boundary

use crate::simulation_parsing::*;
use crate::vessels::*;
use itertools_num::linspace;
use serde::{Deserialize, Serialize};
use splines::{Interpolation, Key, Spline};
use std::cmp::Ordering;
use std::fs::File;
use tracing::{event, info, instrument, span, warn, Level};

/// Contains the simulation parameters and the vessels.
#[derive(Debug, Serialize, Deserialize)]
pub struct Inlet {
    time: Vec<f64>,
    flow: Vec<f64>,
    cardiac_T: f64,
    mid_idx: usize,
    max_idx: usize,
}

impl Inlet {
    pub fn new(data: &Vec<InletRaw>) -> Inlet {
        Inlet {
            time: data.iter().map(|(x, _)| *x).collect(),
            flow: data.iter().map(|(_, y)| *y).collect(),
            cardiac_T: data.last().unwrap().0,
            mid_idx: data.len() / 2,
            max_idx: data.len() - 2,
        }
    }

    pub fn flow_at_time(&self, t: f64) -> f64 {
        let cardiac_nb = (t / self.cardiac_T).trunc();
        let t = t - cardiac_nb * self.cardiac_T;

        let mut step = self.mid_idx;
        let mut idx = self.mid_idx;

        while step > 1 {
            step = (1 + step) / 2;
            idx = match (self.time[idx]).partial_cmp(&t) {
                Some(Ordering::Less) => self.max_idx.min(idx + step),
                Some(Ordering::Greater) => {
                    if idx > step {
                        idx - step
                    } else {
                        0
                    }
                }
                _ => idx,
            };
        }

        if idx < self.max_idx && (self.time[idx] - t).abs() > (self.time[idx + 1] - t).abs() {
            return self.interpolate(idx + 1, t);
        } else if idx > 0 && (self.time[idx] - t).abs() > (self.time[idx - 1] - t).abs() {
            return self.interpolate(idx - 1, t);
        }

        self.interpolate(idx, t)
    }

    fn interpolate(&self, idx: usize, t: f64) -> f64 {
        self.flow[idx]
            + (t - self.time[idx]) * (self.flow[idx + 1] - self.flow[idx]) / (self.time[idx + 1] - self.time[idx])
    }
}
