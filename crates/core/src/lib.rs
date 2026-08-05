//! Shared types and utilities used by the Cotton2K domain crates.

#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(dead_code)]
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
pub mod general_functions;
pub mod global_defs;
pub use global_defs::*;
pub mod input_functions;
pub mod model_state;
pub use model_state::*;
pub mod profile;
pub mod root_impedance;
pub use profile::*;
pub use root_impedance::RootImpedanceTables;
pub mod utils;

#[derive(Debug)]
/// Error produced while preparing or advancing a Cotton2K model state.
///
/// `level == 0` denotes a recoverable simulation stop; higher levels are
/// propagated as failures by the root runner.
pub struct Cotton2KError {
    /// Severity/termination level used by the legacy simulation protocol.
    pub level: u8,
    /// Human-readable failure or stop reason.
    pub message: String,
}

impl std::fmt::Display for Cotton2KError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for Cotton2KError {}

impl Cotton2KError {
    /// Returns the legacy severity level.
    pub fn level(&self) -> u8 {
        self.level
    }

    /// Returns the caller-visible error message.
    pub fn message(&self) -> &str {
        &self.message
    }
}
