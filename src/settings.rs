use structopt::StructOpt;

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
}
