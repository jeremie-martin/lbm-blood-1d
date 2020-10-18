//! Initialize a terminal and a log file [Subscribers](tracing_subscriber::fmt::Subscriber) with a custom [FormatTime](tracing_subscriber::fmt::time::FormatTime).
//!
//! Both the `stdio` and log file [Subscribers](tracing_subscriber::fmt::Subscriber) will write the
//! - Absolute time, formatted with `"%Y-%m-%d %H:%M:%S"`
//! - Elapsed number of seconds since the beginning of the execution

use crate::settings::*;
use std::fmt;
use std::time::Instant;
use tracing::*;

#[derive(Debug, Clone, Eq, PartialEq)]
struct TimeLog {
    epoch: Instant,
    format: String,
}

impl Default for TimeLog {
    fn default() -> Self {
        TimeLog {
            epoch: Instant::now(),
            format: "%Y-%m-%d %H:%M:%S".to_string(),
        }
    }
}

impl tracing_subscriber::fmt::time::FormatTime for TimeLog {
    fn format_time(&self, w: &mut dyn fmt::Write) -> fmt::Result {
        let e = self.epoch.elapsed();
        let time = chrono::Local::now();
        write!(
            w,
            "{} (+{:03}.{:06}s)",
            time.format(&self.format),
            e.as_secs(),
            e.subsec_micros(),
        )
    }
}

/// Installs [subscribers](tracing_subscriber::fmt::Subscriber) given the command line arguments
///
/// Returns the [guards](tracing_appender::non_blocking::WorkerGuard) that we create to keep them alive,
/// which is needed to record events outside of this function
///
/// # Examples
///
/// Basic usage:
///
/// ```
/// tracing_init(Settings { tracing_filter: "info", logs_path: "logs/" })
/// ```
pub fn tracing_init(args: &Settings) -> Vec<impl Drop> {
    use std::fs;
    use tracing_subscriber::prelude::*;
    use tracing_subscriber::{registry, EnvFilter};

    // To hold the guards that we create, they will cause the logs to be flushed when they're dropped.
    let mut _guards = vec![];

    let filter = EnvFilter::try_new(&args.tracing_filter)
        .or_else(|_| EnvFilter::try_from_default_env())
        .or_else(|_| EnvFilter::try_new("info"))
        .unwrap();

    // Create the terminal writer layer.
    let (non_blocking, _stdio_guard) = tracing_appender::non_blocking(std::io::stdout());
    _guards.push(_stdio_guard);

    // Try to create the log file's parent folders.
    let log_folders_created = fs::create_dir_all(&args.logs_path);
    const LOG_FILENAME: &str = "simulation.log";

    let format = tracing_subscriber::fmt::format()
        .with_target(false)
        .with_thread_ids(false)
        .with_thread_names(false)
        .with_timer(TimeLog::default());

    let stdout_layer = tracing_subscriber::fmt::layer()
        .with_writer(non_blocking)
        .event_format(format.clone());

    match log_folders_created {
        // Terminal + file
        Ok(_) => {
            let file_appender = tracing_appender::rolling::never(&args.logs_path, LOG_FILENAME);
            let (non_blocking_file, _file_guard) = tracing_appender::non_blocking(file_appender);
            let file_layer = tracing_subscriber::fmt::layer()
                .with_writer(non_blocking_file)
                .event_format(format.clone());

            _guards.push(_file_guard);
            registry()
                .with(stdout_layer)
                .with(file_layer)
                .with(filter)
                .init();
        }
        // Terminal
        Err(e) => {
            error!(?e, "Failed to create log file!.",);
            registry().with(stdout_layer).with(filter).init();
        }
    };

    // Return the guards
    _guards
}
