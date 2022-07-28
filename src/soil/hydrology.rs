use crate::utils::{fmax, fmin};
use crate::{
    alpha, beta, dl, psiq, qpsi, thad, thetar, thts, wk, ActualTranspiration, AppliedWater,
    AverageSoilPsi, ClayVolumeFraction, DayStart, Daynum, ElCondSatSoilToday, GetFromClim,
    GetTargetStress, Irrig, IrrigMethod, LastIrrigation, LightIntercept, NitrogenUptake,
    NumIrrigations, NumLayersWithRoots, PredictDripIrrigation, PredictSurfaceIrrigation,
    PsiOnTranspiration, PsiOsmotic, ReferenceTransp, RootColNumLeft, RootColNumRight,
    RootWtCapblUptake, RowSpace, SandVolumeFraction, SoilHorizonNum, SoilPsi, SupplyNH4N,
    SupplyNO3N, TotalRequiredN, VolWaterContent, CLIMATE_METRIC_RAIN,
};
use ndarray::prelude::*;
use ndarray::Array;

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
/// Computes the amount of water (mm) applied by a predicted irrigation.
///
/// It is called from SoilProcedures().
///
/// It calls GetTargetStress(), PredictDripIrrigation(), PredictSurfaceIrrigation(),
///
/// The following global variables are referenced here:
///
/// * [AppliedWater]
/// * [Daynum]
/// * [IrrigMethod]
///
/// The following global variable is set here:
/// * [LastIrrigation]
pub unsafe fn ComputeIrrigation() {
    let TargetStress = GetTargetStress();
    if TargetStress == -9999. {
        return;
    }
    if IrrigMethod == 2 {
        PredictDripIrrigation(TargetStress);
    } else {
        PredictSurfaceIrrigation(TargetStress);
    }
    // If the amount of water to be applied (AppliedWater) is non zero update the date of last irrigation, and write report in output file *.B01.
    if AppliedWater > 1e-5 {
        LastIrrigation = Daynum;
    }
}

/// Computes the uptake of water by plant roots from the soil
/// (i.e., actual transpiration rate). It is called from SoilProcedures().
///
/// It calls PsiOnTranspiration(), psiq(), PsiOsmotic().
///
/// The following global variables are referenced:
///
/// dl, LightIntercept, nk, nl, NumLayersWithRoots, ReferenceTransp,
/// RootColNumLeft, RootColNumRight, RowSpace, SoilHorizonNum, thetar,
/// TotalRequiredN, wk.
///
/// The following global variables are set:
/// ActualTranspiration, SoilPsi, VolWaterContent.
pub unsafe fn WaterUptake() {
    // Compute the modified light interception factor (LightInter1) for use in computing transpiration rate.
    // modified light interception factor by canopy
    let LightInter1 = fmin(fmax(LightIntercept * 1.55 - 0.32, LightIntercept), 1.);
    // The potential transpiration is the product of the daytime Penman equation and LightInter1.
    let PotentialTranspiration = ReferenceTransp * LightInter1;
    // uptake factor, computed as a ratio, for each soil cell
    let mut upf = Array::zeros([40, 20]);
    // actual transpiration from each soil cell, cm3 per day
    let mut uptk = Array::<f64, _>::zeros([40, 20]);
    // sum of actual transpiration from all soil soil cells, cm3 per day.
    let mut sumep = 0.;
    // Compute the reduction due to soil moisture supply by function PsiOnTranspiration().
    // the actual transpiration converted to cm3 per slab units.
    let mut Transp = 0.10 * RowSpace * PotentialTranspiration * PsiOnTranspiration(AverageSoilPsi);
    // the cumulative difference between computed transpiration and actual transpiration, in cm3, due to limitation of PWP.
    let mut difupt;

    loop {
        let mut supf = 0.; // sum of upf for all soil cells
        for l in 0..NumLayersWithRoots as usize {
            let j = SoilHorizonNum[l] as usize;
            // Compute, for each layer, the lower and upper water content limits for the transpiration function.
            // These are set from limiting soil water potentials (-15 to -1 bars).
            // lower limit of water content for the transpiration function
            let vh2lo;
            // upper limit of water content for the transpiration function
            let vh2hi;
            vh2lo = qpsi(-15., thad[l], thts[l], alpha[j], beta[j]);
            vh2hi = qpsi(-1., thad[l], thts[l], alpha[j], beta[j]);
            for k in RootColNumLeft[l] as usize..RootColNumRight[l] as usize + 1 {
                // reduction factor for water uptake, caused by low levels of soil water, as a linear function of VolWaterContent, between vh2lo and vh2hi.
                let redfac = fmin(
                    fmax((VolWaterContent[l][k] - vh2lo) / (vh2hi - vh2lo), 0.),
                    1.,
                );
                // The computed 'uptake factor' (upf) for each soil cell is the product of 'root weight capable of uptake' and redfac.
                upf.slice_mut(s![l, k])
                    .fill(RootWtCapblUptake[l][k] * redfac);
                supf += upf.slice(s![l, k]).first().unwrap();
            }
        }
        difupt = 0.;
        for l in 0..NumLayersWithRoots as usize {
            for k in RootColNumLeft[l] as usize..RootColNumRight[l] as usize + 1 {
                if upf.slice(s![l, k]).first().unwrap() > &0. && VolWaterContent[l][k] > thetar[l] {
                    // The amount of water extracted from each cell is proportional to its 'uptake factor'.
                    // transpiration from a soil cell, cm3 per day
                    let mut upth2o = Transp * upf.slice(s![l, k]).first().unwrap() / supf;
                    // Update VolWaterContent cell, storing its previous value as vh2ocx.
                    // previous value of VolWaterContent of this cell
                    let vh2ocx = VolWaterContent[l][k];
                    VolWaterContent[l][k] -= upth2o / (dl[l] * wk[k]);
                    // If the new value of VolWaterContent is less than the permanent wilting point, modify the value of upth2o so that VolWaterContent will be equal to it.
                    if VolWaterContent[l][k] < thetar[l] {
                        VolWaterContent[l][k] = thetar[l];
                        // Compute the difference due to this correction and add it to difupt.
                        // intermediate computation of upth2o
                        let xupt = (vh2ocx - thetar[l]) * dl[l] * wk[k];
                        difupt += upth2o - xupt;
                        upth2o = xupt;
                    }
                    if upth2o < 0. {
                        upth2o = 0.;
                    }
                    // Compute sumep as the sum of the actual amount of water extracted from all soil cells.
                    // Recalculate uptk of this soil cell as cumulative upth2o.
                    sumep += upth2o;
                    uptk.slice_mut(s![l, k]).mapv_inplace(|x| x + upth2o);
                }
            }
        }
        // If difupt is greater than zero, redefine the variable Transp as difuptfor use in next loop.
        if difupt > 0. {
            Transp = difupt;
        } else {
            break;
        }
    }
    // recompute SoilPsi for all soil cells with roots by calling function PSIQ,
    for l in 0..NumLayersWithRoots as usize {
        let j = SoilHorizonNum[l] as usize;
        for k in RootColNumLeft[l] as usize..RootColNumRight[l] as usize + 1 {
            SoilPsi[l][k] = psiq(VolWaterContent[l][k], thad[l], thts[l], alpha[j], beta[j])
                - PsiOsmotic(VolWaterContent[l][k], thts[l], ElCondSatSoilToday);
        }
    }
    // compute ActualTranspiration as actual water transpired, in mm.
    ActualTranspiration = sumep * 10. / RowSpace;
    // Zeroize the amounts of NH4 and NO3 nitrogen taken up from the soil.
    SupplyNH4N = 0.;
    SupplyNO3N = 0.;
    // Compute the proportional N requirement from each soil cell with roots, and call function NitrogenUptake() to compute nitrogen uptake.
    if sumep > 0. && TotalRequiredN > 0. {
        for l in 0..NumLayersWithRoots as usize {
            for k in RootColNumLeft[l] as usize..RootColNumRight[l] as usize + 1 {
                if uptk.slice(s![l, k]).first().unwrap() > &0. {
                    // proportional allocation of TotalRequiredN to each cell
                    let reqnc = TotalRequiredN * uptk.slice(s![l, k]).first().unwrap() / sumep;
                    NitrogenUptake(l as i32, k as i32, reqnc);
                }
            }
        }
    }
}
