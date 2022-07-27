pub mod hydrology;
mod thermology;

use crate::soil::hydrology::SoilHydrology;
use crate::soil::thermology::SoilThermology;

#[derive(Debug, Clone, Copy)]
pub struct Soil {
    pub thermology: SoilThermology,
    pub hydrology: SoilHydrology,
}

impl Soil {
    pub fn new() -> Self {
        Soil {
            thermology: SoilThermology::new(),
            hydrology: SoilHydrology::new(),
        }
    }
}
