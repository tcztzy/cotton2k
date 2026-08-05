//! Soil sub-model for Cotton2K.

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
#![allow(clippy::new_without_default)]

pub use cotton2k_atmosphere as atmosphere;
pub use cotton2k_core::*;
pub use cotton2k_core::{general_functions, input_functions, model_state, profile, utils};

pub mod soil;
pub use soil::*;
