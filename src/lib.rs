#![warn(missing_docs)]
#![warn(missing_doc_code_examples)]

pub mod logging;
pub mod settings;
pub mod simulation;
pub mod vessels;

use color_eyre::eyre::{self, WrapErr};
use tracing::{event, info, instrument, span, warn, Level};

#[cfg(test)]
mod tests {
    use indoc::indoc;
    use pretty_assertions::{assert_eq, assert_ne};

    use super::*;
}
