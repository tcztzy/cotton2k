use crate::{
    ClayVolumeFraction, DayStart, Daynum, GetFromClim, Irrig, NumIrrigations, SandVolumeFraction,
    CLIMATE_METRIC_RAIN,
};

#[derive(Debug, Clone, Copy)]
enum RunoffPotential {
    Low,
    Moderate,
    High,
}

#[derive(Debug, Clone, Copy)]
pub struct SoilHydrology {
    runoff_potential: RunoffPotential,
    pub runoff: f64,
}

impl SoilHydrology {
    pub fn new() -> Self {
        // The following is computed only the first time the function is called.
        // Infiltration rate is estimated from the percent sand and percent clay in the Ap layer.
        // If clay content is greater than 35%, the soil is assumed to have a higher runoff potential,
        // if clay content is less than 15% and sand is greater than 70%, a lower runoff potential is assumed.
        // Other soils (loams) assumed moderate runoff potential.
        // No 'impermeable' (group D) soils are assumed.
        // References: Schwab, Brady.
        let runoff_potential =
            if unsafe { SandVolumeFraction[0] > 0.70 && ClayVolumeFraction[0] < 0.15 } {
                // Soil group A = 1, low runoff potential
                RunoffPotential::Low
            } else if unsafe { ClayVolumeFraction[0] > 0.35 } {
                // Soil group C = 3, high runoff potential
                RunoffPotential::High
            } else {
                // Soil group B = 2, moderate runoff potential
                RunoffPotential::Moderate
            };
        SoilHydrology {
            runoff_potential,
            runoff: 0.,
        }
    }
    /// This function is called from DayClim() and is executed on each day with raifall more than 2 mm.
    /// It computes the runoff and the retained portion of the rainfall.
    ///
    /// Note: This function is based on the code of GOSSYM. No changes have been made from the original GOSSYM code (except translation to C++).
    /// It has not been validated by actual field measurement.
    ///
    /// It calculates the portion of rainfall that is lost to runoff, and reduces rainfall to the amount which is actually infiltrated into the soil.
    /// It uses the soil conservation service method of estimating runoff.
    /// References:
    /// - Brady, Nyle C. 1984. The nature and properties of soils, 9th ed. Macmillan Publishing Co.
    /// - Schwab, Frevert, Edminster, and Barnes. 1981. Soil and water conservation engineering, 3rd ed. John Wiley & Sons, Inc.
    ///
    /// The following global variables are referenced here:
    /// ClayVolumeFraction, Daynum, DayStart, Irrig (structure), NumIrrigations, SandVolumeFraction.
    ///
    /// The argument used here:  rain = today,s rainfall.
    /// The return value:  the amount of water (mm) lost by runoff.
    pub fn runoff(&self, rain: f64) -> f64 {
        // Adjustment of curve number for soil groups A,B,C.
        let d01 = match self.runoff_potential {
            RunoffPotential::Low => 1.0,
            RunoffPotential::Moderate => 1.09,
            RunoffPotential::High => 1.14,
        };
        // Loop to accumulate 5-day antecedent rainfall (mm) which will affect the soil's ability to accept new rainfall. This also includes all irrigations.
        let mut i01 = unsafe { Daynum - 5 };
        if i01 < unsafe { DayStart } {
            i01 = unsafe { DayStart };
        }
        let i02 = unsafe { Daynum };
        let mut previous_wetting = 0.; // five day total (before this day) of rain and irrigation, mm
        for Dayn in i01..i02 {
            let mut amtirr = 0.; // mm water applied on this day by irrigation
            for i in 0..unsafe { NumIrrigations } as usize {
                if Dayn == unsafe { Irrig[i].day } {
                    amtirr = unsafe { Irrig[i].amount };
                }
            }
            previous_wetting += amtirr + unsafe { GetFromClim(CLIMATE_METRIC_RAIN, Dayn) };
        }
        // Adjusting curve number for antecedent rainfall conditions.
        let d02: f64 = if previous_wetting < 36. {
            // low moisture, low runoff potential.
            match self.runoff_potential {
                RunoffPotential::Low => 0.71,
                RunoffPotential::Moderate => 0.78,
                RunoffPotential::High => 0.83,
            }
        } else if previous_wetting > 53. {
            // wet conditions, high runoff potential.
            match self.runoff_potential {
                RunoffPotential::Low => 1.24,
                RunoffPotential::Moderate => 1.15,
                RunoffPotential::High => 1.1,
            }
        } else {
            // moderate conditions
            1.
        };
        // Assuming straight rows, and good cropping practice:
        // Runoff curve number, adjusted for moisture and soil type.
        let crvnum = 78.0 * d01 * d02;
        // maximum potential difference between rainfall and runoff.
        let d03 = 25400. / crvnum - 254.;
        if rain <= 0.2 * d03 {
            0.
        } else {
            (rain - 0.2 * d03).powi(2) / (rain + 0.8 * d03)
        }
    }
}
