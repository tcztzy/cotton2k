//! Plant sub-model composition and per-plant state.
//!
//! [`Plant`] groups residual carbon, nitrogen allocation, and height used by
//! the growth, phenology, root, and abscission modules. [`Plant::new`] creates
//! the model's initial plant state; daily routines in child modules exchange
//! detailed values through the synchronized legacy state mirror.

pub mod abscission;
pub mod growth;
mod nitrogen;
pub mod phenology;
pub mod root;

use nitrogen::PlantNitrogen;

#[derive(Debug, Clone, Copy)]
/// Compact per-plant scratch state used by the daily growth engine.
pub struct Plant {
    /// residual available carbon for root growth from previous day.
    pub pavail: f64,
    /// Nitrogen demand, reserve, and uptake values for the current day.
    pub nitrogen: PlantNitrogen,
    /// Current canopy height in model distance units.
    pub height: f64,
}

impl Plant {
    /// Creates the model's initial plant scratch state.
    pub fn new() -> Self {
        Plant {
            pavail: 0.,
            nitrogen: PlantNitrogen::new(),
            height: 4.,
        }
    }
}

/// Resets plant growth and phenology scratch state before a new run.
pub fn reset_scratch_state() {
    growth::reset_scratch_state();
    phenology::reset_scratch_state();
}
