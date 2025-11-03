pub mod hydrology;
mod thermology;

use crate::profile::Profile;
use crate::soil::hydrology::SoilHydrology;
use crate::soil::thermology::SoilThermology;

#[derive(Debug, Clone, Copy)]
pub struct Soil {
    pub thermology: SoilThermology,
    pub hydrology: SoilHydrology,
}

impl Soil {
    pub fn new(profile: &Profile) -> Self {
        Soil {
            thermology: SoilThermology::new(profile),
            hydrology: SoilHydrology::new(),
        }
    }
}
