//! # Cotton2K
//!
//! Simulation model for cotton

#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(dead_code)]
// Temporary migration boundary for translated legacy model code. Keep new Rust
// entry points and tests idiomatic, and shrink this list as modules are ported.
#![allow(unused_imports)]
#![allow(unused_assignments)]
#![allow(clippy::approx_constant)]
#![allow(clippy::assign_op_pattern)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::collapsible_match)]
#![allow(clippy::derivable_impls)]
#![allow(clippy::doc_lazy_continuation)]
#![allow(clippy::empty_line_after_doc_comments)]
#![allow(clippy::field_reassign_with_default)]
#![allow(clippy::implicit_saturating_sub)]
#![allow(clippy::items_after_test_module)]
#![allow(clippy::manual_clamp)]
#![allow(clippy::manual_is_multiple_of)]
#![allow(clippy::manual_memcpy)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::needless_return)]
#![allow(clippy::single_match)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::type_complexity)]
#![allow(clippy::unnecessary_cast)]
mod bindings;
pub use bindings::*;
mod atmosphere;
mod general_functions;
mod global_defs;
pub use global_defs::*;
mod input_functions;
mod model_state;
pub use model_state::*;
mod plant;
mod profile;
mod runner;
mod soil;
mod state;
mod utils;
use chrono::NaiveDate;

pub use profile::{Profile, SoilHydraulic, WeatherRecord};
pub use runner::{
    run_job, RunError, RunErrorCode, RunEvent, RunFailure, RunRequest, RunStatus, RunSummary,
};
pub use state::State;

#[derive(Debug)]
pub struct Cotton2KError {
    level: u8,
    message: String,
}

impl std::fmt::Display for Cotton2KError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for Cotton2KError {}

impl Cotton2KError {
    pub fn level(&self) -> u8 {
        self.level
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}
