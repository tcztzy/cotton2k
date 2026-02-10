use crate::general_functions::{psiq, qpsi, wcond, GetFromClim, PsiOsmotic};
use crate::utils::{fmax, fmin};
use crate::{
    airdr, alpha, beta, conmax, dl, nk, nl, noitr, thad, thetar, thetas, thts, wk,
    ActualSoilEvaporation, ActualTranspiration, AppliedWater, AverageSoilPsi, ClayVolumeFraction,
    CumWaterDrained, DayStart, DayStartPredIrrig, DayStopPredIrrig, Daynum, ElCondSatSoilToday,
    FirstBloom, FirstSquare, Irrig, IrrigMethod, IrrigationDepth, Kday, LastIrrigation,
    LightIntercept, MaxIrrigation, MaxWaterCapacity, MinDaysBetweenIrrig, NO3FlowFraction,
    NumGreenBolls, NumIrrigations, NumLayersWithRoots, NumOpenBolls, PerPlantArea, PoreSpace,
    RatioImplicit, ReferenceTransp, RootColNumLeft, RootColNumRight, RootWtCapblUptake, RowSpace,
    SandVolumeFraction, SaturatedHydCond, SoilHorizonNum, SoilNitrogenLoss, SoilPsi, SupplyNH4N,
    SupplyNO3N, TotalRequiredN, TotalSoilNh4N, TotalSoilNitrogen, TotalSoilNo3N, TotalSoilUreaN,
    TotalSoilWater, VolNh4NContent, VolNo3NContent, VolUreaNContent, VolWaterContent, WaterStress,
    WaterTableLayer, CLIMATE_METRIC_RAIN,
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
    /// Note: This function is based on the code of GOSSYM and keeps the original behavior after the Rust rewrite.
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
        let irrig = Irrig.read().expect("Irrig lock poisoned");
        for Dayn in i01..i02 {
            let mut amtirr = 0.; // mm water applied on this day by irrigation
            for i in 0..unsafe { NumIrrigations } as usize {
                if Dayn == irrig[i].day {
                    amtirr = irrig[i].amount;
                }
            }
            previous_wetting += amtirr + GetFromClim(CLIMATE_METRIC_RAIN, Dayn);
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

static mut IRR1ST: bool = false;
static mut REQUIRED_WATER: f64 = 0.0;
static mut N_DAYS_BELOW_TARGET_STRESS: i32 = 0;
static mut N_IRR_LAYERS: i32 = 0;
static mut CAPILLARY_NUMITER: i64 = 0;

unsafe fn get_target_stress() -> f64 {
    const STRESS_TARGET: [f64; 10] = [0.70, 0.95, 0.99, 0.99, 0.99, 0.95, 0.90, 0.80, 0.60, 0.40];
    let mut target_stress;
    if Kday > 0 && FirstSquare <= 0 {
        target_stress = STRESS_TARGET[0];
    } else if FirstBloom <= 0 {
        target_stress = STRESS_TARGET[1];
    } else if Daynum <= FirstBloom + 20 {
        target_stress = STRESS_TARGET[2];
    } else if Daynum <= FirstBloom + 40 {
        target_stress = STRESS_TARGET[3];
    } else if NumOpenBolls <= 0.01 {
        target_stress = STRESS_TARGET[4];
    } else if NumOpenBolls < 0.25 * NumGreenBolls {
        target_stress = STRESS_TARGET[5];
    } else if NumOpenBolls < 0.667 * NumGreenBolls {
        target_stress = STRESS_TARGET[6];
    } else if NumOpenBolls < 1.5 * NumGreenBolls {
        target_stress = STRESS_TARGET[7];
    } else if NumOpenBolls < 4.0 * NumGreenBolls {
        target_stress = STRESS_TARGET[8];
    } else if NumOpenBolls < 9.0 * NumGreenBolls {
        target_stress = STRESS_TARGET[9];
    } else {
        DayStopPredIrrig = Daynum;
        target_stress = -9999.0;
    }

    if target_stress <= 0.0 {
        DayStopPredIrrig = Daynum;
        target_stress = -9999.0;
    }
    target_stress
}

unsafe fn predict_drip_irrigation(target_stress: f64) {
    if Daynum <= DayStartPredIrrig {
        IRR1ST = false;
    }

    if !IRR1ST {
        let mut is_irrigated_today = false;
        let irrig = Irrig.read().expect("Irrig lock poisoned");
        for irrigation_index in 0..NumIrrigations as usize {
            if irrig[irrigation_index].day == Daynum
                || GetFromClim(CLIMATE_METRIC_RAIN, Daynum) > 1.0
            {
                is_irrigated_today = true;
                break;
            }
        }

        if !is_irrigated_today && WaterStress <= 0.99 {
            AppliedWater = fmin(30.0, MaxIrrigation);
            IRR1ST = true;
            REQUIRED_WATER = 0.0;
        }
        return;
    }

    REQUIRED_WATER +=
        ActualTranspiration + ActualSoilEvaporation - GetFromClim(CLIMATE_METRIC_RAIN, Daynum);
    if REQUIRED_WATER < 0.0 {
        REQUIRED_WATER = 0.0;
    }

    if (Daynum - MinDaysBetweenIrrig) >= LastIrrigation {
        let mut irrigation_factor = if target_stress > WaterStress {
            1.20 * target_stress / WaterStress
        } else {
            0.90 * target_stress / WaterStress
        };
        irrigation_factor = irrigation_factor.clamp(0.80, 1.25);

        if REQUIRED_WATER * irrigation_factor > MaxIrrigation {
            AppliedWater = MaxIrrigation;
            REQUIRED_WATER -= MaxIrrigation;
        } else {
            AppliedWater = REQUIRED_WATER * irrigation_factor;
            REQUIRED_WATER = 0.0;
        }
    }
}

unsafe fn predict_surface_irrigation(target_stress: f64) {
    if Daynum <= DayStartPredIrrig {
        N_DAYS_BELOW_TARGET_STRESS = 0;
        let mut accumulated_depth = 0.0;
        for layer in 0..nl as usize {
            accumulated_depth += dl[layer];
            if accumulated_depth > IrrigationDepth {
                N_IRR_LAYERS = layer as i32;
                break;
            }
        }
    }

    if (Daynum - MinDaysBetweenIrrig) >= (LastIrrigation - 2)
        && Daynum > DayStartPredIrrig
        && WaterStress < target_stress
    {
        N_DAYS_BELOW_TARGET_STRESS += 1;
        if N_DAYS_BELOW_TARGET_STRESS >= 3 {
            let mut required_water = 0.0;
            for layer in 0..N_IRR_LAYERS as usize {
                for column in 0..nk as usize {
                    let deficit = MaxWaterCapacity[layer] - VolWaterContent[layer][column];
                    required_water += dl[layer] * wk[column] * deficit;
                }
            }
            AppliedWater = required_water * 10.0 / RowSpace;
            if AppliedWater > MaxIrrigation {
                AppliedWater = MaxIrrigation;
            }
            N_DAYS_BELOW_TARGET_STRESS = 0;
        }
    }
}

pub unsafe fn average_psi() -> f64 {
    const VRCUMIN: f64 = 0.1e-9;
    const VRCUMAX: f64 = 0.025;

    let mut psinum = [0.0; 9];
    let mut sumwat = [0.0; 9];
    let mut sumdl = [0.0; 9];

    for layer in 0..NumLayersWithRoots as usize {
        let horizon = SoilHorizonNum[layer] as usize;
        sumdl[horizon] += dl[layer];
        for column in RootColNumLeft[layer] as usize..=RootColNumRight[layer] as usize {
            if RootWtCapblUptake[layer][column] >= VRCUMIN {
                let weight =
                    dl[layer] * wk[column] * fmin(RootWtCapblUptake[layer][column], VRCUMAX);
                sumwat[horizon] += VolWaterContent[layer][column] * weight;
                psinum[horizon] += weight;
            }
        }
    }

    let mut sumpsi = 0.0;
    let mut sumnum = 0.0;
    for horizon in 0..9 {
        if psinum[horizon] > 0.0 && sumdl[horizon] > 0.0 {
            let avgwat = sumwat[horizon] / psinum[horizon];
            let avgpsi = psiq(
                avgwat,
                airdr[horizon],
                thetas[horizon],
                alpha[horizon],
                beta[horizon],
            ) - PsiOsmotic(avgwat, thetas[horizon], ElCondSatSoilToday);
            sumpsi += avgpsi * psinum[horizon];
            sumnum += psinum[horizon];
        }
    }

    if sumnum > 0.0 {
        sumpsi / sumnum
    } else {
        0.0
    }
}

unsafe fn psi_on_transpiration(psi_average: f64) -> f64 {
    const A: f64 = 20.0;
    const B: f64 = 14.0;
    const C: f64 = 1.0;
    const D: f64 = 0.05;

    let mut effect = ((A + psi_average) / B).powi(3);
    if effect > C {
        effect = C;
    }
    if effect < D {
        effect = D;
    }
    effect
}

unsafe fn nitrogen_uptake(layer: i32, column: i32, reqnc: f64) {
    const HALFN: f64 = 0.08;
    const CPARUPMAX: f64 = 0.5;
    const P1: f64 = 100.0;
    const P2: f64 = 5.0;

    let l = layer as usize;
    let k = column as usize;
    let coeff = 10.0 * RowSpace / (PerPlantArea * dl[l] * wk[k]);

    if VolNo3NContent[l][k] > 0.0 {
        let mut uptake_no3 =
            reqnc * VolNo3NContent[l][k] / (HALFN * VolWaterContent[l][k] + VolNo3NContent[l][k]);
        let uptake_max = CPARUPMAX * VolNo3NContent[l][k];
        if coeff * uptake_no3 < uptake_max {
            VolNo3NContent[l][k] -= coeff * uptake_no3;
        } else {
            VolNo3NContent[l][k] -= uptake_max;
            uptake_no3 = uptake_max / coeff;
        }
        SupplyNO3N += uptake_no3;
    }

    if VolNh4NContent[l][k] > 0.0 {
        let bb = P1 + P2 * VolWaterContent[l][k] - VolNh4NContent[l][k];
        let cc = P2 * VolWaterContent[l][k] * VolNh4NContent[l][k];
        let mut ee = bb * bb + 4.0 * cc;
        if ee < 0.0 {
            ee = 0.0;
        }

        let ammonium_dissolved = 0.5 * (ee.sqrt() - bb);
        if ammonium_dissolved > 0.0 {
            let mut uptake_nh4 =
                reqnc * ammonium_dissolved / (HALFN * VolWaterContent[l][k] + ammonium_dissolved);
            let uptake_max = CPARUPMAX * VolNh4NContent[l][k];
            if coeff * uptake_nh4 < uptake_max {
                VolNh4NContent[l][k] -= coeff * uptake_nh4;
            } else {
                VolNh4NContent[l][k] -= uptake_max;
                uptake_nh4 = uptake_max / coeff;
            }
            SupplyNH4N += uptake_nh4;
        }
    }
}

pub unsafe fn soil_sum() {
    TotalSoilWater = 0.0;
    TotalSoilNo3N = 0.0;
    TotalSoilNh4N = 0.0;
    TotalSoilUreaN = 0.0;

    for layer in 0..nl as usize {
        for column in 0..nk as usize {
            TotalSoilNo3N += VolNo3NContent[layer][column] * dl[layer] * wk[column];
            TotalSoilNh4N += VolNh4NContent[layer][column] * dl[layer] * wk[column];
            TotalSoilUreaN += VolUreaNContent[layer][column] * dl[layer] * wk[column];
            TotalSoilWater += VolWaterContent[layer][column] * dl[layer] * wk[column];
        }
    }

    TotalSoilNitrogen = TotalSoilNo3N + TotalSoilNh4N + TotalSoilUreaN;
    TotalSoilWater = TotalSoilWater * 10.0 / RowSpace;
}

unsafe fn water_balance(q1: &mut [f64; 40], qx: &[f64; 40], dd: &[f64; 40], nn: usize) {
    let mut dev = 0.0;
    let mut dabs = 0.0;
    for i in 0..nn {
        dev += dd[i] * (q1[i] - qx[i]);
        dabs += (q1[i] - qx[i]).abs();
    }
    if dabs > 0.0 {
        for i in 0..nn {
            q1[i] -= (q1[i] - qx[i]).abs() * dev / (dabs * dd[i]);
        }
    }
}

unsafe fn nitrogen_flow(
    nn: usize,
    q01: &[f64; 40],
    q1: &[f64; 40],
    dd: &[f64; 40],
    nit: &mut [f64; 40],
    nur: &mut [f64; 40],
) {
    for i in 0..nn {
        if nur[i] < 1e-20 {
            nur[i] = 0.0;
        }
        if nit[i] < 1e-20 {
            nit[i] = 0.0;
        }
    }

    let mut qdn = [0.0; 40];
    let mut qup = [0.0; 40];
    let mut udn = [0.0; 40];
    let mut uup = [0.0; 40];

    for i in 0..nn {
        let aq0 = q01[i] * dd[i];
        let aq1 = q1[i] * dd[i];
        if i == 0 {
            qup[i] = 0.0;
            uup[i] = 0.0;
        } else {
            qup[i] = -qdn[i - 1];
            uup[i] = -udn[i - 1];
        }

        if i == nn - 1 {
            qdn[i] = 0.0;
            udn[i] = 0.0;
        } else {
            qdn[i] = (aq1 - aq0) * nit[i + 1] / q01[i + 1];
            if qdn[i] < (-0.2 * nit[i] * dd[i]) {
                qdn[i] = -0.2 * nit[i] * dd[i];
            }
            if qdn[i] > (0.2 * nit[i + 1] * dd[i + 1]) {
                qdn[i] = 0.2 * nit[i + 1] * dd[i + 1];
            }

            udn[i] = (aq1 - aq0) * nur[i + 1] / q01[i + 1];
            if udn[i] < (-0.2 * nur[i] * dd[i]) {
                udn[i] = -0.2 * nur[i] * dd[i];
            }
            if udn[i] > (0.2 * nur[i + 1] * dd[i + 1]) {
                udn[i] = 0.2 * nur[i + 1] * dd[i + 1];
            }
        }
    }

    for i in 0..nn {
        nit[i] += (qdn[i] + qup[i]) / dd[i];
        nur[i] += (udn[i] + uup[i]) / dd[i];
    }
}

unsafe fn water_flux(
    q1: &mut [f64; 40],
    psi1: &mut [f64; 40],
    dd: &[f64; 40],
    qr1: &[f64; 40],
    qs1: &[f64; 40],
    pp1: &[f64; 40],
    nn: usize,
    iv: i32,
    ll: usize,
    numiter: i64,
) {
    if nn < 2 {
        return;
    }

    let delt = 1.0 / noitr as f64;
    let mut cond = [0.0; 40];
    let mut kx = [0.0; 40];
    let mut ky = [0.0; 40];

    let mut j = SoilHorizonNum[ll] as usize;
    for i in 0..nn {
        if iv == 1 {
            j = SoilHorizonNum[i] as usize;
        }
        cond[i] = wcond(q1[i], qr1[i], qs1[i], beta[j], SaturatedHydCond[j], pp1[i]);
        kx[i] = 0.0;
        ky[i] = 0.0;
    }

    let mut dy = [0.0; 40];
    let mut avcond = [0.0; 40];
    const COND_MIN: f64 = 0.000006;
    for i in 1..nn {
        dy[i] = 0.5 * (dd[i - 1] + dd[i]);
        if cond[i - 1] <= COND_MIN && cond[i] <= COND_MIN {
            avcond[i] = COND_MIN;
        } else if cond[i - 1] <= COND_MIN && cond[i] > COND_MIN {
            avcond[i] = 2.0 * COND_MIN * cond[i] / (COND_MIN + cond[i]);
        } else if cond[i] <= COND_MIN && cond[i - 1] > COND_MIN {
            avcond[i] = 2.0 * COND_MIN * cond[i - 1] / (COND_MIN + cond[i - 1]);
        } else {
            avcond[i] = 2.0 * cond[i - 1] * cond[i] / (cond[i - 1] + cond[i]);
        }
    }

    let mut qx = [0.0; 40];
    let mut addq = [0.0; 40];
    for i in 0..nn {
        qx[i] = q1[i];
    }

    for i in 1..(nn - 1) {
        let mut deltpsi = (psi1[i - 1] - psi1[i]).clamp(-1000.0, 1000.0);
        if iv == 1 {
            deltpsi += 0.001 * dy[i];
        }
        let mut dumm1 = 1000.0 * avcond[i] * delt / dy[i];
        if dumm1 > conmax * dy[i] {
            dumm1 = conmax * dy[i];
        }

        let mut dqq1 = (1.0 - RatioImplicit) * deltpsi * dumm1;
        let mut deltq = qx[i - 1] - qx[i];
        if dqq1.abs() > (0.25 * deltq).abs() {
            if (deltq > 0.0 && dqq1 < 0.0) || (deltq < 0.0 && dqq1 > 0.0) {
                dqq1 = 0.0;
            } else {
                dqq1 = 0.25 * deltq;
            }
        }

        deltpsi = (psi1[i + 1] - psi1[i]).clamp(-1000.0, 1000.0);
        deltq = qx[i + 1] - qx[i];
        if iv == 1 {
            deltpsi -= 0.001 * dy[i + 1];
        }
        dumm1 = 1000.0 * avcond[i + 1] * delt / dy[i + 1];
        if dumm1 > conmax * dy[i + 1] {
            dumm1 = conmax * dy[i + 1];
        }

        let mut dqq2 = (1.0 - RatioImplicit) * deltpsi * dumm1;
        if dqq2.abs() > (0.25 * deltq).abs() {
            if (deltq > 0.0 && dqq2 < 0.0) || (deltq < 0.0 && dqq2 > 0.0) {
                dqq2 = 0.0;
            } else {
                dqq2 = 0.25 * deltq;
            }
        }

        addq[i] = (dqq1 + dqq2) / dd[i];
        if i == 1 {
            addq[0] = -dqq1 / dd[0];
        }
        if i == nn - 2 {
            addq[nn - 1] = -dqq2 / dd[nn - 1];
        }
    }

    for i in 0..nn {
        q1[i] = qx[i] + addq[i];
        if iv == 1 {
            j = SoilHorizonNum[i] as usize;
        }
        psi1[i] = psiq(q1[i], qr1[i], qs1[i], alpha[j], beta[j]);
    }

    for i in 1..nn {
        ky[i] = 1000.0 * avcond[i] * delt / (dy[i] * dd[i]);
        if ky[i] < 0.0000001 {
            ky[i] = 0.0;
        }
        if ky[i] > conmax {
            ky[i] = conmax;
        }
    }
    for i in 0..(nn - 1) {
        kx[i] = 1000.0 * avcond[i + 1] * delt / (dy[i + 1] * dd[i]);
        if kx[i] < 0.0000001 {
            kx[i] = 0.0;
        }
        if kx[i] > conmax {
            kx[i] = conmax;
        }
    }

    let mut a1 = [0.0; 40];
    let mut b1 = [0.0; 40];
    let mut cau = [0.0; 40];
    let mut cc1 = [0.0; 40];
    let mut d1 = [0.0; 40];
    let mut dau = [0.0; 40];

    for i in 0..nn {
        a1[i] = -kx[i] * RatioImplicit;
        b1[i] = 1.0 + RatioImplicit * (kx[i] + ky[i]);
        cc1[i] = -ky[i] * RatioImplicit;
        if iv == 1 {
            j = SoilHorizonNum[i] as usize;
            a1[i] -= 0.001 * kx[i] * RatioImplicit;
            cc1[i] += 0.001 * ky[i] * RatioImplicit;
        }
        d1[i] = psiq(q1[i], qr1[i], qs1[i], alpha[j], beta[j]);
    }

    if numiter % 2 == 0 {
        cau[nn - 1] = psi1[nn - 1];
        dau[nn - 1] = 0.0;
        for i in (1..(nn - 1)).rev() {
            let p = a1[i] * dau[i + 1] + b1[i];
            dau[i] = -cc1[i] / p;
            cau[i] = (d1[i] - a1[i] * cau[i + 1]) / p;
        }
        if iv == 1 {
            j = SoilHorizonNum[0] as usize;
        }
        psi1[0] = psiq(q1[0], qr1[0], qs1[0], alpha[j], beta[j]);
        for i in 1..(nn - 1) {
            if iv == 1 {
                j = SoilHorizonNum[i] as usize;
            }
            psi1[i] = dau[i] * psi1[i - 1] + cau[i];
            q1[i] = qpsi(psi1[i], qr1[i], qs1[i], alpha[j], beta[j]);
        }
    } else {
        cau[0] = psi1[0];
        dau[0] = 0.0;
        for i in 1..(nn - 1) {
            let p = a1[i] * dau[i - 1] + b1[i];
            dau[i] = -cc1[i] / p;
            cau[i] = (d1[i] - a1[i] * cau[i - 1]) / p;
        }
        if iv == 1 {
            j = SoilHorizonNum[nn - 1] as usize;
        }
        psi1[nn - 1] = psiq(q1[nn - 1], qr1[nn - 1], qs1[nn - 1], alpha[j], beta[j]);
        for i in (1..(nn - 1)).rev() {
            if iv == 1 {
                j = SoilHorizonNum[i] as usize;
            }
            psi1[i] = dau[i] * psi1[i + 1] + cau[i];
            q1[i] = qpsi(psi1[i], qr1[i], qs1[i], alpha[j], beta[j]);
        }
    }

    for i in 0..nn {
        q1[i] = q1[i].clamp(qr1[i], qs1[i]);
        if q1[i] > pp1[i] {
            q1[i] = pp1[i];
        }
    }
    water_balance(q1, &qx, dd, nn);
}

pub unsafe fn drain() -> f64 {
    let mut nlx = nl;
    if WaterTableLayer < nlx {
        nlx = WaterTableLayer;
    }
    if nlx <= 0 {
        return 0.0;
    }

    let mut oldvh2oc = [0.0; 20];
    for k in 0..nk as usize {
        oldvh2oc[k] = VolWaterContent[nlx as usize - 1][k];
    }

    for l in 0..(nlx as usize).saturating_sub(1) {
        let mut avwl = 0.0;
        for k in 0..nk as usize {
            avwl += VolWaterContent[l][k] * wk[k] / RowSpace;
            oldvh2oc[k] = VolWaterContent[l][k];
        }

        let uplimit = MaxWaterCapacity[l];
        if avwl > uplimit {
            let mut wmov = avwl - uplimit;
            wmov = wmov * dl[l] / dl[l + 1];
            for k in 0..nk as usize {
                VolWaterContent[l][k] = uplimit;
                VolWaterContent[l + 1][k] += wmov * wk[k] * nk as f64 / RowSpace;

                let qvout = oldvh2oc[k] - uplimit;
                if qvout > 0.0 {
                    let mut nitconc = VolNo3NContent[l][k] / oldvh2oc[k];
                    if nitconc < 1e-30 {
                        nitconc = 0.0;
                    }
                    let mut nurconc = VolUreaNContent[l][k] / oldvh2oc[k];
                    if nurconc < 1e-30 {
                        nurconc = 0.0;
                    }
                    VolNo3NContent[l][k] = VolWaterContent[l][k] * nitconc;
                    VolUreaNContent[l][k] = VolWaterContent[l][k] * nurconc;

                    let vno3mov = qvout * nitconc;
                    VolNo3NContent[l + 1][k] += NO3FlowFraction[l] * vno3mov * dl[l] / dl[l + 1];
                    VolNo3NContent[l][k] += (1.0 - NO3FlowFraction[l]) * vno3mov;

                    let vnurmov = qvout * nurconc;
                    VolUreaNContent[l + 1][k] += NO3FlowFraction[l] * vnurmov * dl[l] / dl[l + 1];
                    VolUreaNContent[l][k] += (1.0 - NO3FlowFraction[l]) * vnurmov;
                }
            }
        } else {
            for k in 0..nk as usize {
                if VolWaterContent[l][k] > uplimit {
                    let wmov = VolWaterContent[l][k] - uplimit;
                    VolWaterContent[l][k] = uplimit;
                    VolWaterContent[l + 1][k] += wmov * dl[l] / dl[l + 1];

                    let mut nitconc = VolNo3NContent[l][k] / oldvh2oc[k];
                    if nitconc < 1e-30 {
                        nitconc = 0.0;
                    }
                    let mut nurconc = VolUreaNContent[l][k] / oldvh2oc[k];
                    if nurconc < 1e-30 {
                        nurconc = 0.0;
                    }
                    VolNo3NContent[l][k] = VolWaterContent[l][k] * nitconc;
                    VolUreaNContent[l][k] = VolWaterContent[l][k] * nurconc;

                    VolNo3NContent[l + 1][k] +=
                        NO3FlowFraction[l] * wmov * nitconc * dl[l] / dl[l + 1];
                    VolUreaNContent[l + 1][k] +=
                        NO3FlowFraction[l] * wmov * nurconc * dl[l] / dl[l + 1];
                    VolNo3NContent[l][k] += (1.0 - NO3FlowFraction[l]) * wmov * nitconc;
                    VolUreaNContent[l][k] += (1.0 - NO3FlowFraction[l]) * wmov * nurconc;
                }
            }
        }
    }

    let mut drainage = 0.0;
    let bottom_layer = nlx as usize - 1;
    for k in 0..nk as usize {
        if VolWaterContent[bottom_layer][k] > MaxWaterCapacity[bottom_layer] {
            drainage += (VolWaterContent[bottom_layer][k] - MaxWaterCapacity[bottom_layer])
                * dl[bottom_layer]
                * wk[k];

            let mut nitconc = VolNo3NContent[bottom_layer][k] / oldvh2oc[k];
            if nitconc < 1e-30 {
                nitconc = 0.0;
            }
            let mut nurconc = VolUreaNContent[bottom_layer][k] / oldvh2oc[k];
            if nurconc < 1e-30 {
                nurconc = 0.0;
            }

            let saven = (VolNo3NContent[bottom_layer][k] + VolUreaNContent[bottom_layer][k])
                * dl[bottom_layer]
                * wk[k];
            VolWaterContent[bottom_layer][k] = MaxWaterCapacity[bottom_layer];
            VolNo3NContent[bottom_layer][k] = nitconc * MaxWaterCapacity[bottom_layer];
            VolUreaNContent[bottom_layer][k] = nurconc * MaxWaterCapacity[bottom_layer];
            SoilNitrogenLoss += saven
                - (VolNo3NContent[bottom_layer][k] + VolUreaNContent[bottom_layer][k])
                    * dl[bottom_layer]
                    * wk[k];
        }
    }

    drainage
}

pub unsafe fn capillary_flow() {
    let mut wk1 = [0.0; 40];
    if Daynum <= DayStart {
        CAPILLARY_NUMITER = 0;
        for layer in 0..nl as usize {
            wk1[layer] = 0.0;
        }
    }

    CAPILLARY_NUMITER += 1;

    for layer in 0..nl as usize {
        let horizon = SoilHorizonNum[layer] as usize;
        for column in 0..nk as usize {
            SoilPsi[layer][column] = psiq(
                VolWaterContent[layer][column],
                thad[layer],
                thts[layer],
                alpha[horizon],
                beta[horizon],
            ) - PsiOsmotic(
                VolWaterContent[layer][column],
                thts[layer],
                ElCondSatSoilToday,
            );
        }
    }

    let mut nlx = nl;
    if WaterTableLayer < nlx {
        nlx = WaterTableLayer - 1;
    }

    let mut q01 = [0.0; 40];
    let mut q1 = [0.0; 40];
    let mut psi1 = [0.0; 40];
    let mut nit = [0.0; 40];
    let mut nur = [0.0; 40];

    let mut dl1 = [0.0; 40];
    let mut thad1 = [0.0; 40];
    let mut thts1 = [0.0; 40];
    let mut ps1 = [0.0; 40];
    for layer in 0..nl as usize {
        dl1[layer] = dl[layer];
        thad1[layer] = thad[layer];
        thts1[layer] = thts[layer];
        ps1[layer] = PoreSpace[layer];
    }

    if nlx > 0 {
        for column in 0..nk as usize {
            for layer in 0..nlx as usize {
                q1[layer] = VolWaterContent[layer][column];
                q01[layer] = VolWaterContent[layer][column];
                psi1[layer] = SoilPsi[layer][column]
                    + PsiOsmotic(
                        VolWaterContent[layer][column],
                        thts[layer],
                        ElCondSatSoilToday,
                    );
                nit[layer] = VolNo3NContent[layer][column];
                nur[layer] = VolUreaNContent[layer][column];
            }

            water_flux(
                &mut q1,
                &mut psi1,
                &dl1,
                &thad1,
                &thts1,
                &ps1,
                nlx as usize,
                1,
                0,
                CAPILLARY_NUMITER,
            );
            nitrogen_flow(nl as usize, &q01, &q1, &dl1, &mut nit, &mut nur);

            for layer in 0..nlx as usize {
                VolWaterContent[layer][column] = q1[layer];
                VolNo3NContent[layer][column] = nit[layer];
                VolUreaNContent[layer][column] = nur[layer];
                SoilPsi[layer][column] = psi1[layer]
                    - PsiOsmotic(
                        VolWaterContent[layer][column],
                        thts[layer],
                        ElCondSatSoilToday,
                    );
            }
        }

        let mut pp1 = [0.0; 40];
        let mut qr1 = [0.0; 40];
        let mut qs1 = [0.0; 40];

        for layer in 0..nlx as usize {
            for column in 0..nk as usize {
                q1[column] = VolWaterContent[layer][column];
                q01[column] = VolWaterContent[layer][column];
                psi1[column] = SoilPsi[layer][column]
                    + PsiOsmotic(
                        VolWaterContent[layer][column],
                        thts[layer],
                        ElCondSatSoilToday,
                    );
                qr1[column] = thad[layer];
                qs1[column] = thts[layer];
                pp1[column] = PoreSpace[layer];
                nit[column] = VolNo3NContent[layer][column];
                nur[column] = VolUreaNContent[layer][column];
                wk1[column] = wk[column];
            }

            water_flux(
                &mut q1,
                &mut psi1,
                &wk1,
                &qr1,
                &qs1,
                &pp1,
                nk as usize,
                0,
                layer,
                CAPILLARY_NUMITER,
            );
            nitrogen_flow(nk as usize, &q01, &q1, &wk1, &mut nit, &mut nur);

            for column in 0..nk as usize {
                VolWaterContent[layer][column] = q1[column];
                SoilPsi[layer][column] = psi1[column]
                    - PsiOsmotic(
                        VolWaterContent[layer][column],
                        thts[layer],
                        ElCondSatSoilToday,
                    );
                VolNo3NContent[layer][column] = nit[column];
                VolUreaNContent[layer][column] = nur[column];
            }
        }
    }

    let water_drained_out = drain();
    if water_drained_out > 0.0 {
        CumWaterDrained += 10.0 * water_drained_out / RowSpace;
    }

    for layer in 0..nl as usize {
        let horizon = SoilHorizonNum[layer] as usize;
        for column in 0..nk as usize {
            SoilPsi[layer][column] = psiq(
                VolWaterContent[layer][column],
                thad[layer],
                thts[layer],
                alpha[horizon],
                beta[horizon],
            ) - PsiOsmotic(
                VolWaterContent[layer][column],
                thts[layer],
                ElCondSatSoilToday,
            );
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
    let target_stress = get_target_stress();
    if target_stress == -9999. {
        return;
    }
    if IrrigMethod == 2 {
        predict_drip_irrigation(target_stress);
    } else {
        predict_surface_irrigation(target_stress);
    }
    // If the amount of water to be applied (AppliedWater) is non zero update the date of last irrigation, and write report in output file *.B01.
    if AppliedWater > 1e-5 {
        LastIrrigation = Daynum;
    }
}

/// Computes the uptake of water by plant roots from the soil
/// (i.e., actual transpiration rate). It is called from SoilProcedures().
///
/// It calls psi_on_transpiration(), psiq(), PsiOsmotic().
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
    // Compute the reduction due to soil moisture supply by function psi_on_transpiration().
    // the actual transpiration converted to cm3 per slab units.
    let mut Transp =
        0.10 * RowSpace * PotentialTranspiration * psi_on_transpiration(AverageSoilPsi);
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
    // Compute the proportional N requirement from each soil cell with roots, and call function nitrogen_uptake() to compute nitrogen uptake.
    if sumep > 0. && TotalRequiredN > 0. {
        for l in 0..NumLayersWithRoots as usize {
            for k in RootColNumLeft[l] as usize..RootColNumRight[l] as usize + 1 {
                if uptk.slice(s![l, k]).first().unwrap() > &0. {
                    // proportional allocation of TotalRequiredN to each cell
                    let reqnc = TotalRequiredN * uptk.slice(s![l, k]).first().unwrap() / sumep;
                    nitrogen_uptake(l as i32, k as i32, reqnc);
                }
            }
        }
    }
}
