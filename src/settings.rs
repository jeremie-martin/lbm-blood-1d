//! Defines the command line arguments

use structopt::StructOpt;

/// Structure used by (structopt)[structopt] to parse the command line arguments
#[derive(Debug, StructOpt)]
pub struct Settings {
    /// Tracing filter.
    ///
    /// Can be any of "error", "warn", "info", "debug", or
    /// "trace". Supports more granular filtering, as well; see documentation for
    /// [`tracing_subscriber::EnvFilter`][EnvFilter].
    ///
    /// [EnvFilter]: https://docs.rs/tracing-subscriber/latest/tracing_subscriber/struct.EnvFilter.html
    #[structopt(long, default_value = "info")]
    pub tracing_filter: String,

    /// Logs directory
    #[structopt(long, default_value = "logs")]
    pub logs_path: String,

    /// Save every `x` second
    #[structopt(short = "x", long, default_value = "0.01")]
    pub time_between_save: f64,

    /// Space step [m]
    #[structopt(long, default_value = "0.01")]
    pub dx: f64,

    #[structopt(long, default_value = "16.666")]
    pub fact: f64,
}
