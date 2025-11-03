pub mod hydrology;
mod thermodynamics;

use crate::profile::Profile;
use crate::soil::hydrology::SoilHydrology;
use crate::soil::thermodynamics::SoilThermodynamics;

#[derive(Debug, Clone, Copy)]
pub struct Soil {
    pub thermodynamics: SoilThermodynamics,
    pub hydrology: SoilHydrology,
}

impl Soil {
    pub fn new(profile: &Profile) -> Self {
        Soil {
            thermodynamics: SoilThermodynamics::new(profile),
            hydrology: SoilHydrology::new(),
        }
    }
}
