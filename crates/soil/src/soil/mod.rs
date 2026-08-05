//! Soil sub-model composition and daily scratch-state reset boundary.
//!
//! [`Soil`] groups hydrology and thermodynamics for a simulation state. The
//! child modules add nitrogen and temperature transitions; [`Soil::new`]
//! initializes their per-run configuration from [`Profile`]. All daily soil
//! routines exchange data through the shared legacy state and perform no I/O.

pub mod hydrology;
pub mod nitrogen;
mod temperature_abi;
mod thermodynamics;

use crate::profile::Profile;
use crate::soil::hydrology::SoilHydrology;
use crate::soil::thermodynamics::SoilThermodynamics;

#[derive(Debug, Clone, Copy)]
/// Per-day soil-domain state combining thermal and hydrological scratch data.
pub struct Soil {
    /// Layer temperature and surface-energy-balance scratch state.
    pub thermodynamics: SoilThermodynamics,
    /// Water movement, runoff, drainage, and uptake scratch state.
    pub hydrology: SoilHydrology,
}

impl Soil {
    /// Initializes soil sub-model state from an initialized profile.
    pub fn new(profile: &Profile) -> Self {
        Soil {
            thermodynamics: SoilThermodynamics::new(profile),
            hydrology: SoilHydrology::new(),
        }
    }
}

/// Resets process-wide scratch state before starting a new simulation.
pub fn reset_scratch_state() {
    hydrology::reset_scratch_state();
    nitrogen::reset_scratch_state();
    temperature_abi::reset_scratch_state();
}
