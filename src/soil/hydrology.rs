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
use ndarray::Array2;

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

fn get_target_stress() -> f64 {
    const STRESS_TARGET: [f64; 10] = [0.70, 0.95, 0.99, 0.99, 0.99, 0.95, 0.90, 0.80, 0.60, 0.40];

    let (kday, first_square, first_bloom, daynum, num_open_bolls, num_green_bolls) = unsafe {
        (
            Kday,
            FirstSquare,
            FirstBloom,
            Daynum,
            NumOpenBolls,
            NumGreenBolls,
        )
    };

    let mut stop_prediction = false;
    let mut target_stress;
    if kday > 0 && first_square <= 0 {
        target_stress = STRESS_TARGET[0];
    } else if first_bloom <= 0 {
        target_stress = STRESS_TARGET[1];
    } else if daynum <= first_bloom + 20 {
        target_stress = STRESS_TARGET[2];
    } else if daynum <= first_bloom + 40 {
        target_stress = STRESS_TARGET[3];
    } else if num_open_bolls <= 0.01 {
        target_stress = STRESS_TARGET[4];
    } else if num_open_bolls < 0.25 * num_green_bolls {
        target_stress = STRESS_TARGET[5];
    } else if num_open_bolls < 0.667 * num_green_bolls {
        target_stress = STRESS_TARGET[6];
    } else if num_open_bolls < 1.5 * num_green_bolls {
        target_stress = STRESS_TARGET[7];
    } else if num_open_bolls < 4.0 * num_green_bolls {
        target_stress = STRESS_TARGET[8];
    } else if num_open_bolls < 9.0 * num_green_bolls {
        target_stress = STRESS_TARGET[9];
    } else {
        stop_prediction = true;
        target_stress = -9999.0;
    }

    if target_stress <= 0.0 {
        stop_prediction = true;
        target_stress = -9999.0;
    }

    if stop_prediction {
        unsafe {
            DayStopPredIrrig = daynum;
        }
    }

    target_stress
}

fn predict_drip_irrigation(target_stress: f64) {
    let (
        daynum,
        day_start_pred_irrig,
        num_irrigations,
        water_stress,
        max_irrigation,
        min_days_between_irrig,
        last_irrigation,
        actual_transpiration,
        actual_soil_evaporation,
        mut irr1st,
        mut required_water,
    ) = unsafe {
        (
            Daynum,
            DayStartPredIrrig,
            NumIrrigations as usize,
            WaterStress,
            MaxIrrigation,
            MinDaysBetweenIrrig,
            LastIrrigation,
            ActualTranspiration,
            ActualSoilEvaporation,
            IRR1ST,
            REQUIRED_WATER,
        )
    };
    let mut applied_water = None;

    if daynum <= day_start_pred_irrig {
        irr1st = false;
    }

    if !irr1st {
        let mut is_irrigated_today = false;
        let irrig = Irrig.read().expect("Irrig lock poisoned");
        for irrigation_index in 0..num_irrigations {
            if irrig[irrigation_index].day == daynum
                || GetFromClim(CLIMATE_METRIC_RAIN, daynum) > 1.0
            {
                is_irrigated_today = true;
                break;
            }
        }

        if !is_irrigated_today && water_stress <= 0.99 {
            applied_water = Some(fmin(30.0, max_irrigation));
            irr1st = true;
            required_water = 0.0;
        }
    } else {
        required_water += actual_transpiration + actual_soil_evaporation
            - GetFromClim(CLIMATE_METRIC_RAIN, daynum);
        if required_water < 0.0 {
            required_water = 0.0;
        }

        if (daynum - min_days_between_irrig) >= last_irrigation {
            let mut irrigation_factor = if target_stress > water_stress {
                1.20 * target_stress / water_stress
            } else {
                0.90 * target_stress / water_stress
            };
            irrigation_factor = irrigation_factor.clamp(0.80, 1.25);

            if required_water * irrigation_factor > max_irrigation {
                applied_water = Some(max_irrigation);
                required_water -= max_irrigation;
            } else {
                applied_water = Some(required_water * irrigation_factor);
                required_water = 0.0;
            }
        }
    }

    unsafe {
        IRR1ST = irr1st;
        REQUIRED_WATER = required_water;
        if let Some(amount) = applied_water {
            AppliedWater = amount;
        }
    }
}

fn predict_surface_irrigation(target_stress: f64) {
    let (
        daynum,
        day_start_pred_irrig,
        min_days_between_irrig,
        last_irrigation,
        water_stress,
        irrigation_depth,
        max_irrigation,
        row_space,
        num_cols,
        dl_layer,
        wk_col,
        max_water_capacity,
        vol_water,
        mut n_days_below_target_stress,
        mut n_irr_layers,
    ) = unsafe {
        let num_layers = nl as usize;
        let num_cols = nk as usize;
        let mut dl_layer = vec![0.0f64; num_layers];
        let mut max_water_capacity = vec![0.0f64; num_layers];
        let mut wk_col = vec![0.0f64; num_cols];
        let mut vol_water = Array2::<f64>::zeros((num_layers, num_cols));

        for layer in 0..num_layers {
            dl_layer[layer] = dl[layer];
            max_water_capacity[layer] = MaxWaterCapacity[layer];
            for column in 0..num_cols {
                vol_water[[layer, column]] = VolWaterContent[layer][column];
            }
        }
        for column in 0..num_cols {
            wk_col[column] = wk[column];
        }

        (
            Daynum,
            DayStartPredIrrig,
            MinDaysBetweenIrrig,
            LastIrrigation,
            WaterStress,
            IrrigationDepth,
            MaxIrrigation,
            RowSpace,
            num_cols,
            dl_layer,
            wk_col,
            max_water_capacity,
            vol_water,
            N_DAYS_BELOW_TARGET_STRESS,
            N_IRR_LAYERS,
        )
    };
    let mut applied_water = None;

    if daynum <= day_start_pred_irrig {
        n_days_below_target_stress = 0;
        let mut accumulated_depth = 0.0;
        for (layer, depth) in dl_layer.iter().enumerate() {
            accumulated_depth += depth;
            if accumulated_depth > irrigation_depth {
                n_irr_layers = layer as i32;
                break;
            }
        }
    }

    if (daynum - min_days_between_irrig) >= (last_irrigation - 2)
        && daynum > day_start_pred_irrig
        && water_stress < target_stress
    {
        n_days_below_target_stress += 1;
        if n_days_below_target_stress >= 3 {
            let mut required_water = 0.0;
            for layer in 0..n_irr_layers.max(0) as usize {
                for column in 0..num_cols {
                    let deficit = max_water_capacity[layer] - vol_water[[layer, column]];
                    required_water += dl_layer[layer] * wk_col[column] * deficit;
                }
            }

            let mut amount = required_water * 10.0 / row_space;
            if amount > max_irrigation {
                amount = max_irrigation;
            }
            applied_water = Some(amount);
            n_days_below_target_stress = 0;
        }
    }

    unsafe {
        N_DAYS_BELOW_TARGET_STRESS = n_days_below_target_stress;
        N_IRR_LAYERS = n_irr_layers;
        if let Some(amount) = applied_water {
            AppliedWater = amount;
        }
    }
}

pub fn average_psi() -> f64 {
    const VRCUMIN: f64 = 0.1e-9;
    const VRCUMAX: f64 = 0.025;

    let mut psinum = [0.0; 9];
    let mut sumwat = [0.0; 9];
    let mut sumdl = [0.0; 9];

    let (
        num_layers,
        root_left,
        root_right,
        soil_horizon,
        dl_layer,
        wk_col,
        root_uptake,
        vol_water,
        airdr_horizon,
        thetas_horizon,
        alpha_horizon,
        beta_horizon,
        el_cond_sat_soil,
    ) = unsafe {
        let num_layers = NumLayersWithRoots as usize;
        let mut root_left = vec![0usize; num_layers];
        let mut root_right = vec![0usize; num_layers];
        let mut soil_horizon = vec![0usize; num_layers];
        let mut dl_layer = vec![0.0f64; num_layers];
        for layer in 0..num_layers {
            root_left[layer] = RootColNumLeft[layer] as usize;
            root_right[layer] = RootColNumRight[layer] as usize;
            soil_horizon[layer] = SoilHorizonNum[layer] as usize;
            dl_layer[layer] = dl[layer];
        }

        let num_cols = nk as usize;
        let mut wk_col = vec![0.0f64; num_cols];
        for column in 0..num_cols {
            wk_col[column] = wk[column];
        }

        let mut root_uptake = Array2::<f64>::zeros((num_layers, num_cols));
        let mut vol_water = Array2::<f64>::zeros((num_layers, num_cols));
        for layer in 0..num_layers {
            for column in 0..num_cols {
                root_uptake[[layer, column]] = RootWtCapblUptake[layer][column];
                vol_water[[layer, column]] = VolWaterContent[layer][column];
            }
        }

        let mut airdr_horizon = [0.0f64; 9];
        let mut thetas_horizon = [0.0f64; 9];
        let mut alpha_horizon = [0.0f64; 9];
        let mut beta_horizon = [0.0f64; 9];
        for horizon in 0..9 {
            airdr_horizon[horizon] = airdr[horizon];
            thetas_horizon[horizon] = thetas[horizon];
            alpha_horizon[horizon] = alpha[horizon];
            beta_horizon[horizon] = beta[horizon];
        }
        let el_cond_sat_soil = ElCondSatSoilToday;

        (
            num_layers,
            root_left,
            root_right,
            soil_horizon,
            dl_layer,
            wk_col,
            root_uptake,
            vol_water,
            airdr_horizon,
            thetas_horizon,
            alpha_horizon,
            beta_horizon,
            el_cond_sat_soil,
        )
    };

    for layer in 0..num_layers {
        let horizon = soil_horizon[layer];
        sumdl[horizon] += dl_layer[layer];
        for column in root_left[layer]..=root_right[layer] {
            let uptake = root_uptake[[layer, column]];
            if uptake >= VRCUMIN {
                let weight = dl_layer[layer] * wk_col[column] * fmin(uptake, VRCUMAX);
                sumwat[horizon] += vol_water[[layer, column]] * weight;
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
                airdr_horizon[horizon],
                thetas_horizon[horizon],
                alpha_horizon[horizon],
                beta_horizon[horizon],
            ) - PsiOsmotic(avgwat, thetas_horizon[horizon], el_cond_sat_soil);
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

fn psi_on_transpiration(psi_average: f64) -> f64 {
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

fn nitrogen_uptake(
    layer: usize,
    column: usize,
    reqnc: f64,
    row_space: f64,
    per_plant_area: f64,
    dl_layer: &[f64],
    wk_col: &[f64],
    vol_water: &Array2<f64>,
    vol_no3: &mut Array2<f64>,
    vol_nh4: &mut Array2<f64>,
    supply_no3: &mut f64,
    supply_nh4: &mut f64,
) {
    const HALFN: f64 = 0.08;
    const CPARUPMAX: f64 = 0.5;
    const P1: f64 = 100.0;
    const P2: f64 = 5.0;

    let coeff = 10.0 * row_space / (per_plant_area * dl_layer[layer] * wk_col[column]);
    let water = vol_water[[layer, column]];

    if vol_no3[[layer, column]] > 0.0 {
        let mut uptake_no3 =
            reqnc * vol_no3[[layer, column]] / (HALFN * water + vol_no3[[layer, column]]);
        let uptake_max = CPARUPMAX * vol_no3[[layer, column]];
        if coeff * uptake_no3 < uptake_max {
            vol_no3[[layer, column]] -= coeff * uptake_no3;
        } else {
            vol_no3[[layer, column]] -= uptake_max;
            uptake_no3 = uptake_max / coeff;
        }
        *supply_no3 += uptake_no3;
    }

    if vol_nh4[[layer, column]] > 0.0 {
        let bb = P1 + P2 * water - vol_nh4[[layer, column]];
        let cc = P2 * water * vol_nh4[[layer, column]];
        let mut ee = bb * bb + 4.0 * cc;
        if ee < 0.0 {
            ee = 0.0;
        }

        let ammonium_dissolved = 0.5 * (ee.sqrt() - bb);
        if ammonium_dissolved > 0.0 {
            let mut uptake_nh4 = reqnc * ammonium_dissolved / (HALFN * water + ammonium_dissolved);
            let uptake_max = CPARUPMAX * vol_nh4[[layer, column]];
            if coeff * uptake_nh4 < uptake_max {
                vol_nh4[[layer, column]] -= coeff * uptake_nh4;
            } else {
                vol_nh4[[layer, column]] -= uptake_max;
                uptake_nh4 = uptake_max / coeff;
            }
            *supply_nh4 += uptake_nh4;
        }
    }
}

pub fn soil_sum() {
    let (num_layers, num_cols, row_space, dl_layer, wk_col, vol_no3, vol_nh4, vol_urea, vol_water) = unsafe {
        let num_layers = nl as usize;
        let num_cols = nk as usize;
        let row_space = RowSpace;

        let mut dl_layer = vec![0.0f64; num_layers];
        for layer in 0..num_layers {
            dl_layer[layer] = dl[layer];
        }
        let mut wk_col = vec![0.0f64; num_cols];
        for column in 0..num_cols {
            wk_col[column] = wk[column];
        }

        let mut vol_no3 = Array2::<f64>::zeros((num_layers, num_cols));
        let mut vol_nh4 = Array2::<f64>::zeros((num_layers, num_cols));
        let mut vol_urea = Array2::<f64>::zeros((num_layers, num_cols));
        let mut vol_water = Array2::<f64>::zeros((num_layers, num_cols));
        for layer in 0..num_layers {
            for column in 0..num_cols {
                vol_no3[[layer, column]] = VolNo3NContent[layer][column];
                vol_nh4[[layer, column]] = VolNh4NContent[layer][column];
                vol_urea[[layer, column]] = VolUreaNContent[layer][column];
                vol_water[[layer, column]] = VolWaterContent[layer][column];
            }
        }

        (
            num_layers, num_cols, row_space, dl_layer, wk_col, vol_no3, vol_nh4, vol_urea,
            vol_water,
        )
    };

    let cell_area = Array2::from_shape_fn((num_layers, num_cols), |(layer, column)| {
        dl_layer[layer] * wk_col[column]
    });
    let total_soil_no3 = (&vol_no3 * &cell_area).sum();
    let total_soil_nh4 = (&vol_nh4 * &cell_area).sum();
    let total_soil_urea = (&vol_urea * &cell_area).sum();
    let total_soil_water = (&vol_water * &cell_area).sum();

    unsafe {
        TotalSoilNo3N = total_soil_no3;
        TotalSoilNh4N = total_soil_nh4;
        TotalSoilUreaN = total_soil_urea;
        TotalSoilNitrogen = TotalSoilNo3N + TotalSoilNh4N + TotalSoilUreaN;
        TotalSoilWater = total_soil_water * 10.0 / row_space;
    }
}

fn water_balance(q1: &mut [f64; 40], qx: &[f64; 40], dd: &[f64; 40], nn: usize) {
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

fn nitrogen_flow(
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
pub fn ComputeIrrigation() {
    let target_stress = get_target_stress();
    if target_stress == -9999. {
        return;
    }
    unsafe {
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
pub fn WaterUptake() {
    let (
        num_layers,
        num_cols,
        row_space,
        light_intercept,
        reference_transp,
        average_soil_psi,
        total_required_n,
        per_plant_area,
        el_cond_sat_soil,
        dl_layer,
        wk_col,
        thad_layer,
        thts_layer,
        thetar_layer,
        soil_horizon,
        root_left,
        root_right,
        alpha_horizon,
        beta_horizon,
        root_uptake,
        mut vol_water,
        mut soil_psi,
        mut vol_no3,
        mut vol_nh4,
    ) = unsafe {
        let num_layers = NumLayersWithRoots as usize;
        let num_cols = nk as usize;
        let row_space = RowSpace;
        let light_intercept = LightIntercept;
        let reference_transp = ReferenceTransp;
        let average_soil_psi = AverageSoilPsi;
        let total_required_n = TotalRequiredN;
        let per_plant_area = PerPlantArea;
        let el_cond_sat_soil = ElCondSatSoilToday;

        let mut dl_layer = vec![0.0f64; num_layers];
        let mut thad_layer = vec![0.0f64; num_layers];
        let mut thts_layer = vec![0.0f64; num_layers];
        let mut thetar_layer = vec![0.0f64; num_layers];
        let mut soil_horizon = vec![0usize; num_layers];
        let mut root_left = vec![0usize; num_layers];
        let mut root_right = vec![0usize; num_layers];
        for layer in 0..num_layers {
            dl_layer[layer] = dl[layer];
            thad_layer[layer] = thad[layer];
            thts_layer[layer] = thts[layer];
            thetar_layer[layer] = thetar[layer];
            soil_horizon[layer] = SoilHorizonNum[layer] as usize;
            root_left[layer] = RootColNumLeft[layer] as usize;
            root_right[layer] = RootColNumRight[layer] as usize;
        }

        let mut wk_col = vec![0.0f64; num_cols];
        for column in 0..num_cols {
            wk_col[column] = wk[column];
        }

        let mut alpha_horizon = [0.0f64; 9];
        let mut beta_horizon = [0.0f64; 9];
        for horizon in 0..9 {
            alpha_horizon[horizon] = alpha[horizon];
            beta_horizon[horizon] = beta[horizon];
        }

        let mut root_uptake = Array2::<f64>::zeros((num_layers, num_cols));
        let mut vol_water = Array2::<f64>::zeros((num_layers, num_cols));
        let mut soil_psi = Array2::<f64>::zeros((num_layers, num_cols));
        let mut vol_no3 = Array2::<f64>::zeros((num_layers, num_cols));
        let mut vol_nh4 = Array2::<f64>::zeros((num_layers, num_cols));
        for layer in 0..num_layers {
            for column in 0..num_cols {
                root_uptake[[layer, column]] = RootWtCapblUptake[layer][column];
                vol_water[[layer, column]] = VolWaterContent[layer][column];
                soil_psi[[layer, column]] = SoilPsi[layer][column];
                vol_no3[[layer, column]] = VolNo3NContent[layer][column];
                vol_nh4[[layer, column]] = VolNh4NContent[layer][column];
            }
        }

        (
            num_layers,
            num_cols,
            row_space,
            light_intercept,
            reference_transp,
            average_soil_psi,
            total_required_n,
            per_plant_area,
            el_cond_sat_soil,
            dl_layer,
            wk_col,
            thad_layer,
            thts_layer,
            thetar_layer,
            soil_horizon,
            root_left,
            root_right,
            alpha_horizon,
            beta_horizon,
            root_uptake,
            vol_water,
            soil_psi,
            vol_no3,
            vol_nh4,
        )
    };

    // Compute the modified light interception factor (LightInter1) for use in computing transpiration rate.
    let light_inter1 = fmin(fmax(light_intercept * 1.55 - 0.32, light_intercept), 1.);
    let potential_transpiration = reference_transp * light_inter1;

    // uptake factor, computed as a ratio, for each soil cell
    let mut upf = Array2::<f64>::zeros((num_layers, num_cols));
    // actual transpiration from each soil cell, cm3 per day
    let mut uptk = Array2::<f64>::zeros((num_layers, num_cols));
    // sum of actual transpiration from all soil cells, cm3 per day.
    let mut sumep = 0.;
    let mut transp =
        0.10 * row_space * potential_transpiration * psi_on_transpiration(average_soil_psi);

    loop {
        let mut supf = 0.;
        for layer in 0..num_layers {
            let horizon = soil_horizon[layer];
            let vh2lo = qpsi(
                -15.,
                thad_layer[layer],
                thts_layer[layer],
                alpha_horizon[horizon],
                beta_horizon[horizon],
            );
            let vh2hi = qpsi(
                -1.,
                thad_layer[layer],
                thts_layer[layer],
                alpha_horizon[horizon],
                beta_horizon[horizon],
            );
            for column in root_left[layer]..=root_right[layer] {
                let redfac = fmin(
                    fmax((vol_water[[layer, column]] - vh2lo) / (vh2hi - vh2lo), 0.),
                    1.,
                );
                let uptake = root_uptake[[layer, column]] * redfac;
                upf[[layer, column]] = uptake;
                supf += uptake;
            }
        }

        let mut difupt = 0.;
        for layer in 0..num_layers {
            for column in root_left[layer]..=root_right[layer] {
                let uptake_factor = upf[[layer, column]];
                if uptake_factor > 0. && vol_water[[layer, column]] > thetar_layer[layer] {
                    let mut upth2o = transp * uptake_factor / supf;
                    let vh2ocx = vol_water[[layer, column]];
                    vol_water[[layer, column]] -= upth2o / (dl_layer[layer] * wk_col[column]);
                    if vol_water[[layer, column]] < thetar_layer[layer] {
                        vol_water[[layer, column]] = thetar_layer[layer];
                        let xupt =
                            (vh2ocx - thetar_layer[layer]) * dl_layer[layer] * wk_col[column];
                        difupt += upth2o - xupt;
                        upth2o = xupt;
                    }
                    if upth2o < 0. {
                        upth2o = 0.;
                    }
                    sumep += upth2o;
                    uptk[[layer, column]] += upth2o;
                }
            }
        }

        if difupt > 0. {
            transp = difupt;
        } else {
            break;
        }
    }

    for layer in 0..num_layers {
        let horizon = soil_horizon[layer];
        for column in root_left[layer]..=root_right[layer] {
            soil_psi[[layer, column]] = psiq(
                vol_water[[layer, column]],
                thad_layer[layer],
                thts_layer[layer],
                alpha_horizon[horizon],
                beta_horizon[horizon],
            ) - PsiOsmotic(
                vol_water[[layer, column]],
                thts_layer[layer],
                el_cond_sat_soil,
            );
        }
    }

    let actual_transpiration = sumep * 10. / row_space;
    let mut supply_no3 = 0.0;
    let mut supply_nh4 = 0.0;
    if sumep > 0. && total_required_n > 0. {
        for layer in 0..num_layers {
            for column in root_left[layer]..=root_right[layer] {
                if uptk[[layer, column]] > 0. {
                    let reqnc = total_required_n * uptk[[layer, column]] / sumep;
                    nitrogen_uptake(
                        layer,
                        column,
                        reqnc,
                        row_space,
                        per_plant_area,
                        &dl_layer,
                        &wk_col,
                        &vol_water,
                        &mut vol_no3,
                        &mut vol_nh4,
                        &mut supply_no3,
                        &mut supply_nh4,
                    );
                }
            }
        }
    }

    unsafe {
        for layer in 0..num_layers {
            for column in 0..num_cols {
                VolWaterContent[layer][column] = vol_water[[layer, column]];
                SoilPsi[layer][column] = soil_psi[[layer, column]];
                VolNo3NContent[layer][column] = vol_no3[[layer, column]];
                VolNh4NContent[layer][column] = vol_nh4[[layer, column]];
            }
        }
        ActualTranspiration = actual_transpiration;
        SupplyNO3N = supply_no3;
        SupplyNH4N = supply_nh4;
    }
}
