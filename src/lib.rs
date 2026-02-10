//! # Cotton2K
//!
//! Simulation model for cotton

#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(dead_code)]
mod bindings;
pub use bindings::*;
mod atmosphere;
mod general_functions;
mod global_defs;
pub use global_defs::*;
mod input_functions;
mod plant;
mod profile;
mod soil;
mod state;
mod utils;
use chrono::NaiveDate;

pub use profile::{Profile, SoilHydraulic, WeatherRecord};
pub use state::State;

#[derive(Debug)]
pub struct Cotton2KError {
    level: u8,
    message: String,
}

impl std::fmt::Display for Cotton2KError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", "Cotton2KError")
    }
}

impl std::error::Error for Cotton2KError {}
