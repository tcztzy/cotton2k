pub mod abscission;
pub mod growth;
mod nitrogen;
pub mod phenology;
pub mod root;

use nitrogen::PlantNitrogen;

#[derive(Debug, Clone, Copy)]
pub struct Plant {
    /// residual available carbon for root growth from previous day.
    pub pavail: f64,
    pub nitrogen: PlantNitrogen,
    pub height: f64,
}

impl Plant {
    pub fn new() -> Self {
        Plant {
            pavail: 0.,
            nitrogen: PlantNitrogen::new(),
            height: 4.,
        }
    }
}

pub(crate) fn reset_scratch_state() {
    growth::reset_scratch_state();
    phenology::reset_scratch_state();
}
