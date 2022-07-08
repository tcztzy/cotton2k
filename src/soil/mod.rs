mod thermology;

use crate::soil::thermology::SoilThermology;

#[derive(Debug, Clone, Copy)]
pub struct Soil {
    pub thermology: SoilThermology,
}

impl Soil {
    pub fn new() -> Self {
        Soil {
            thermology: SoilThermology::new(),
        }
    }
}
