#![warn(missing_docs)]
#![warn(missing_doc_code_examples)]
#![allow(non_snake_case)]

//! One-dimensional blood flow dynamics simulation using the Lattice Boltzmann Method

pub mod boundary_conditions;
pub mod constants;
pub mod logging;
pub mod populations;
pub mod settings;
pub mod simulation;
pub mod simulation_parsing;
pub mod vessels;
pub mod vessels_cells;
pub mod vessels_parsing;
