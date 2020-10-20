extern crate chrono;
extern crate csv;
extern crate itertools_num;
extern crate replace_with;
extern crate ryu; // float to String scientific notation
extern crate serde;
extern crate serde_json;
extern crate splines;

use color_eyre::eyre;
use lbm_blood_1d::logging::tracing_init;
use lbm_blood_1d::settings::Settings;
use lbm_blood_1d::simulation::*;
use structopt::StructOpt;
use tracing::*;

#[instrument]
fn main() -> eyre::Result<()> {
    let args = Settings::from_args();
    let _guards = tracing_init(&args);
    color_eyre::install()?;

    info!("Execution started.");

    Simulation::new("vascularNetworks/bifur.json");

    Ok(())
}
