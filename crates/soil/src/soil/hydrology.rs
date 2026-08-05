//! Soil water movement, irrigation, runoff, drainage, and root uptake.
//!
//! [`SoilHydrology`] holds hydrology scratch configuration while the public
//! routines advance water and nitrogen fluxes in the synchronized legacy grid.
//! Inputs are climate, management, and soil state; outputs are updated water
//! contents, potentials, flux totals, and irrigation decisions. The routines
//! are stateful, non-I/O, and must be called in the engine's daily order.

use crate::general_functions::{psiq, qpsi, wcond, GetFromClim, PsiOsmotic};
use crate::model_state::{for_each_layer_col_span, for_each_row};
use crate::utils::{fmax, fmin};
use crate::LegacyGlobalState;
use crate::{
    airdr, alpha, beta, conmax, dl, nk, nl, noitr, thad, thetar, thetas, thts, wk,
    ActualSoilEvaporation, ActualTranspiration, AppliedWater, AverageSoilPsi, CumWaterDrained,
    DayStart, DayStartPredIrrig, DayStopPredIrrig, Daynum, ElCondSatSoilToday, FirstBloom,
    FirstSquare, Irrig, IrrigMethod, IrrigationDepth, Kday, LastIrrigation, LightIntercept,
    MaxIrrigation, MaxWaterCapacity, MinDaysBetweenIrrig, NO3FlowFraction, NumGreenBolls,
    NumIrrigations, NumLayersWithRoots, NumOpenBolls, PerPlantArea, PoreSpace, RatioImplicit,
    ReferenceTransp, RootColNumLeft, RootColNumRight, RootWtCapblUptake, RowSpace,
    SaturatedHydCond, SoilHorizonNum, SoilNitrogenLoss, SoilPsi, SupplyNH4N, SupplyNO3N,
    TotalRequiredN, TotalSoilNh4N, TotalSoilNitrogen, TotalSoilNo3N, TotalSoilUreaN,
    TotalSoilWater, VolNh4NContent, VolNo3NContent, VolUreaNContent, VolWaterContent, WaterStress,
    WaterTableLayer, CLIMATE_METRIC_RAIN,
};
use ndarray::{s, Array1, Array2};
use std::sync::{LazyLock, RwLock};

#[derive(Debug, Clone, Copy)]
enum RunoffPotential {
    Low,
    Moderate,
    High,
}

#[derive(Debug, Clone, Copy)]
/// Per-run hydrology configuration and accumulated runoff.
pub struct SoilHydrology {
    runoff_potential: RunoffPotential,
    /// Rainfall lost as surface runoff during the current calculation.
    pub runoff: f64,
}

impl SoilHydrology {
    /// Determines runoff potential from the initialized surface soil fractions.
    pub fn new() -> Self {
        // The following is computed only the first time the function is called.
        // Infiltration rate is estimated from the percent sand and percent clay in the Ap layer.
        // If clay content is greater than 35%, the soil is assumed to have a higher runoff potential,
        // if clay content is less than 15% and sand is greater than 70%, a lower runoff potential is assumed.
        // Other soils (loams) assumed moderate runoff potential.
        // No 'impermeable' (group D) soils are assumed.
        // References: Schwab, Brady.
        let legacy = LegacyGlobalState::from_globals();
        let surface_sand_fraction = legacy.sand_volume_fraction[0];
        let surface_clay_fraction = legacy.clay_volume_fraction[0];
        let runoff_potential = if surface_sand_fraction > 0.70 && surface_clay_fraction < 0.15 {
            // Soil group A = 1, low runoff potential
            RunoffPotential::Low
        } else if surface_clay_fraction > 0.35 {
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
        let legacy = LegacyGlobalState::from_globals();
        let mut i01 = legacy.daynum - 5;
        let day_start = legacy.day_start;
        let i02 = legacy.daynum;
        let num_irrigations = legacy.num_irrigations as usize;
        if i01 < day_start {
            i01 = day_start;
        }
        let mut previous_wetting = 0.; // five day total (before this day) of rain and irrigation, mm
        let irrig = Irrig.read().expect("Irrig lock poisoned");
        for Dayn in i01..i02 {
            let mut amtirr = 0.; // mm water applied on this day by irrigation
            for i in 0..num_irrigations {
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

#[derive(Debug, Clone, Copy, Default)]
struct HydrologyScratch {
    irr1st: bool,
    required_water: f64,
    n_days_below_target_stress: i32,
    n_irr_layers: i32,
    capillary_numiter: i64,
}

static HYDROLOGY_SCRATCH: LazyLock<RwLock<HydrologyScratch>> =
    LazyLock::new(|| RwLock::new(HydrologyScratch::default()));

fn read_hydrology_scratch() -> HydrologyScratch {
    *HYDROLOGY_SCRATCH
        .read()
        .expect("hydrology scratch state lock should not be poisoned")
}

fn write_hydrology_scratch(scratch: HydrologyScratch) {
    *HYDROLOGY_SCRATCH
        .write()
        .expect("hydrology scratch state lock should not be poisoned") = scratch;
}

pub(crate) fn reset_scratch_state() {
    write_hydrology_scratch(HydrologyScratch::default());
}

fn get_target_stress() -> f64 {
    const STRESS_TARGET: [f64; 10] = [0.70, 0.95, 0.99, 0.99, 0.99, 0.95, 0.90, 0.80, 0.60, 0.40];

    let mut legacy = LegacyGlobalState::from_globals();
    let mut stop_prediction = false;
    let mut target_stress;
    if legacy.kday > 0 && legacy.first_square <= 0 {
        target_stress = STRESS_TARGET[0];
    } else if legacy.first_bloom <= 0 {
        target_stress = STRESS_TARGET[1];
    } else if legacy.daynum <= legacy.first_bloom + 20 {
        target_stress = STRESS_TARGET[2];
    } else if legacy.daynum <= legacy.first_bloom + 40 {
        target_stress = STRESS_TARGET[3];
    } else if legacy.num_open_bolls <= 0.01 {
        target_stress = STRESS_TARGET[4];
    } else if legacy.num_open_bolls < 0.25 * legacy.num_green_bolls {
        target_stress = STRESS_TARGET[5];
    } else if legacy.num_open_bolls < 0.667 * legacy.num_green_bolls {
        target_stress = STRESS_TARGET[6];
    } else if legacy.num_open_bolls < 1.5 * legacy.num_green_bolls {
        target_stress = STRESS_TARGET[7];
    } else if legacy.num_open_bolls < 4.0 * legacy.num_green_bolls {
        target_stress = STRESS_TARGET[8];
    } else if legacy.num_open_bolls < 9.0 * legacy.num_green_bolls {
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
        legacy.day_stop_pred_irrig = legacy.daynum;
        legacy.write_to_globals();
    }
    target_stress
}
fn predict_drip_irrigation(target_stress: f64) {
    let mut legacy = LegacyGlobalState::from_globals();
    let daynum = legacy.daynum;
    let day_start_pred_irrig = legacy.day_start_pred_irrig;
    let num_irrigations = legacy.num_irrigations as usize;
    let water_stress = legacy.water_stress;
    let max_irrigation = legacy.max_irrigation;
    let min_days_between_irrig = legacy.min_days_between_irrig;
    let last_irrigation = legacy.last_irrigation;
    let actual_transpiration = legacy.actual_transpiration;
    let actual_soil_evaporation = legacy.actual_soil_evaporation;
    let mut scratch = read_hydrology_scratch();
    let mut irr1st = scratch.irr1st;
    let mut required_water = scratch.required_water;

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

    scratch.irr1st = irr1st;
    scratch.required_water = required_water;
    write_hydrology_scratch(scratch);
    if let Some(amount) = applied_water {
        legacy.applied_water = amount;
        legacy.write_to_globals();
    }
}
fn predict_surface_irrigation(target_stress: f64) {
    let mut legacy = LegacyGlobalState::from_globals();
    let daynum = legacy.daynum;
    let day_start_pred_irrig = legacy.day_start_pred_irrig;
    let min_days_between_irrig = legacy.min_days_between_irrig;
    let last_irrigation = legacy.last_irrigation;
    let water_stress = legacy.water_stress;
    let irrigation_depth = legacy.irrigation_depth;
    let max_irrigation = legacy.max_irrigation;
    let row_space = legacy.row_space;
    let num_layers = legacy.nl as usize;
    let num_cols = legacy.nk as usize;

    let mut scratch = read_hydrology_scratch();
    let mut n_days_below_target_stress = scratch.n_days_below_target_stress;
    let mut n_irr_layers = scratch.n_irr_layers;
    let mut applied_water = None;

    if daynum <= day_start_pred_irrig {
        n_days_below_target_stress = 0;
        let mut accumulated_depth = 0.0;
        for (layer, depth) in legacy.dl.iter().take(num_layers).enumerate() {
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
            for_each_row(
                &legacy.vol_water_content,
                n_irr_layers.max(0) as usize,
                |layer, row| {
                    for (column, &vol_water) in row.iter().take(num_cols).enumerate() {
                        let deficit = legacy.max_water_capacity[layer] - vol_water;
                        required_water += legacy.dl[layer] * legacy.wk[column] * deficit;
                    }
                },
            );

            let mut amount = required_water * 10.0 / row_space;
            if amount > max_irrigation {
                amount = max_irrigation;
            }
            applied_water = Some(amount);
            n_days_below_target_stress = 0;
        }
    }

    scratch.n_days_below_target_stress = n_days_below_target_stress;
    scratch.n_irr_layers = n_irr_layers;
    write_hydrology_scratch(scratch);
    if let Some(amount) = applied_water {
        legacy.applied_water = amount;
        legacy.write_to_globals();
    }
}
/// Computes the root-zone average soil water potential in model units.
pub fn average_psi() -> f64 {
    const VRCUMIN: f64 = 0.1e-9;
    const VRCUMAX: f64 = 0.025;

    let legacy = LegacyGlobalState::from_globals();
    let mut psinum = [0.0; 9];
    let mut sumwat = [0.0; 9];
    let mut sumdl = [0.0; 9];

    let num_layers = legacy.num_layers_with_roots as usize;

    for layer in 0..num_layers {
        let horizon = legacy.soil_horizon_num[layer] as usize;
        sumdl[horizon] += legacy.dl[layer];
    }
    for_each_layer_col_span(
        num_layers,
        legacy.nk as usize,
        legacy.root_col_num_left.as_slice().unwrap(),
        legacy.root_col_num_right.as_slice().unwrap(),
        |layer, column| {
            let horizon = legacy.soil_horizon_num[layer] as usize;
            let uptake = legacy.root_wt_capbl_uptake[[layer, column]];
            if uptake >= VRCUMIN {
                let weight = legacy.dl[layer] * legacy.wk[column] * fmin(uptake, VRCUMAX);
                sumwat[horizon] += legacy.vol_water_content[[layer, column]] * weight;
                psinum[horizon] += weight;
            }
        },
    );

    let mut sumpsi = 0.0;
    let mut sumnum = 0.0;
    for horizon in 0..9 {
        if psinum[horizon] > 0.0 && sumdl[horizon] > 0.0 {
            let avgwat = sumwat[horizon] / psinum[horizon];
            let avgpsi = psiq(
                avgwat,
                legacy.airdr[horizon],
                legacy.thetas[horizon],
                legacy.alpha[horizon],
                legacy.beta[horizon],
            ) - PsiOsmotic(
                avgwat,
                legacy.thetas[horizon],
                legacy.el_cond_sat_soil_today,
            );
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

/// Recomputes aggregate soil water and nitrogen totals from the grid.
pub fn soil_sum() {
    let mut legacy = LegacyGlobalState::from_globals();
    let num_layers = legacy.nl as usize;
    let num_cols = legacy.nk as usize;

    let mut total_soil_no3 = 0.0;
    let mut total_soil_nh4 = 0.0;
    let mut total_soil_urea = 0.0;
    let mut total_soil_water = 0.0;
    for_each_row(&legacy.vol_water_content, num_layers, |layer, water_row| {
        let no3_row = legacy.vol_no3_n_content.row(layer);
        let nh4_row = legacy.vol_nh4_n_content.row(layer);
        let urea_row = legacy.vol_urea_n_content.row(layer);
        let layer_depth = legacy.dl[layer];
        for column in 0..num_cols {
            let area = layer_depth * legacy.wk[column];
            total_soil_no3 += no3_row[column] * area;
            total_soil_nh4 += nh4_row[column] * area;
            total_soil_urea += urea_row[column] * area;
            total_soil_water += water_row[column] * area;
        }
    });

    legacy.total_soil_no3_n = total_soil_no3;
    legacy.total_soil_nh4_n = total_soil_nh4;
    legacy.total_soil_urea_n = total_soil_urea;
    legacy.total_soil_nitrogen = total_soil_no3 + total_soil_nh4 + total_soil_urea;
    legacy.total_soil_water = total_soil_water * 10.0 / legacy.row_space;
    legacy.write_to_globals();
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

fn water_flux(
    q1: &mut [f64; 40],
    psi1: &mut [f64; 40],
    dd: &[f64; 40],
    qr1: &[f64; 40],
    qs1: &[f64; 40],
    pp1: &[f64; 40],
    nn: usize,
    iv: i32,
    ll: usize,
    noitr_value: i32,
    soil_horizon_num: &Array1<i32>,
    beta_profile: &Array1<f64>,
    saturated_hyd_cond_profile: &Array1<f64>,
    alpha_profile: &Array1<f64>,
    ratio_implicit: f64,
    conmax_value: f64,
    numiter: i64,
) {
    if nn < 2 {
        return;
    }

    let delt = 1.0 / noitr_value as f64;
    let mut cond = [0.0; 40];
    let mut kx = [0.0; 40];
    let mut ky = [0.0; 40];

    let mut j = soil_horizon_num[ll] as usize;
    for i in 0..nn {
        if iv == 1 {
            j = soil_horizon_num[i] as usize;
        }
        cond[i] = wcond(
            q1[i],
            qr1[i],
            qs1[i],
            beta_profile[j],
            saturated_hyd_cond_profile[j],
            pp1[i],
        );
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
        if dumm1 > conmax_value * dy[i] {
            dumm1 = conmax_value * dy[i];
        }

        let mut dqq1 = (1.0 - ratio_implicit) * deltpsi * dumm1;
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
        if dumm1 > conmax_value * dy[i + 1] {
            dumm1 = conmax_value * dy[i + 1];
        }

        let mut dqq2 = (1.0 - ratio_implicit) * deltpsi * dumm1;
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
            j = soil_horizon_num[i] as usize;
        }
        psi1[i] = psiq(q1[i], qr1[i], qs1[i], alpha_profile[j], beta_profile[j]);
    }

    for i in 1..nn {
        ky[i] = 1000.0 * avcond[i] * delt / (dy[i] * dd[i]);
        if ky[i] < 0.0000001 {
            ky[i] = 0.0;
        }
        if ky[i] > conmax_value {
            ky[i] = conmax_value;
        }
    }
    for i in 0..(nn - 1) {
        kx[i] = 1000.0 * avcond[i + 1] * delt / (dy[i + 1] * dd[i]);
        if kx[i] < 0.0000001 {
            kx[i] = 0.0;
        }
        if kx[i] > conmax_value {
            kx[i] = conmax_value;
        }
    }

    let mut a1 = [0.0; 40];
    let mut b1 = [0.0; 40];
    let mut cau = [0.0; 40];
    let mut cc1 = [0.0; 40];
    let mut d1 = [0.0; 40];
    let mut dau = [0.0; 40];

    for i in 0..nn {
        a1[i] = -kx[i] * ratio_implicit;
        b1[i] = 1.0 + ratio_implicit * (kx[i] + ky[i]);
        cc1[i] = -ky[i] * ratio_implicit;
        if iv == 1 {
            j = soil_horizon_num[i] as usize;
            a1[i] -= 0.001 * kx[i] * ratio_implicit;
            cc1[i] += 0.001 * ky[i] * ratio_implicit;
        }
        d1[i] = psiq(q1[i], qr1[i], qs1[i], alpha_profile[j], beta_profile[j]);
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
            j = soil_horizon_num[0] as usize;
        }
        psi1[0] = psiq(q1[0], qr1[0], qs1[0], alpha_profile[j], beta_profile[j]);
        for i in 1..(nn - 1) {
            if iv == 1 {
                j = soil_horizon_num[i] as usize;
            }
            psi1[i] = dau[i] * psi1[i - 1] + cau[i];
            q1[i] = qpsi(psi1[i], qr1[i], qs1[i], alpha_profile[j], beta_profile[j]);
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
            j = soil_horizon_num[nn - 1] as usize;
        }
        psi1[nn - 1] = psiq(
            q1[nn - 1],
            qr1[nn - 1],
            qs1[nn - 1],
            alpha_profile[j],
            beta_profile[j],
        );
        for i in (1..(nn - 1)).rev() {
            if iv == 1 {
                j = soil_horizon_num[i] as usize;
            }
            psi1[i] = dau[i] * psi1[i + 1] + cau[i];
            q1[i] = qpsi(psi1[i], qr1[i], qs1[i], alpha_profile[j], beta_profile[j]);
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

/// Drains excess water below the model's free-drainage limit.
pub fn drain() -> f64 {
    let mut legacy = LegacyGlobalState::from_globals();
    let drainage = drain_with_legacy(&mut legacy);
    legacy.write_to_globals();
    drainage
}

fn drain_with_legacy(legacy: &mut LegacyGlobalState) -> f64 {
    fn nutrient_concentration(content: f64, water_content: f64) -> f64 {
        let concentration = content / water_content;
        if concentration < 1e-30 {
            0.0
        } else {
            concentration
        }
    }

    fn move_mobile_nutrients(
        moved_water: f64,
        nitconc: f64,
        nurconc: f64,
        no3_flow_fraction: f64,
        layer_depth_ratio: f64,
    ) -> ((f64, f64), (f64, f64)) {
        let vno3mov = moved_water * nitconc;
        let vnurmov = moved_water * nurconc;
        (
            (
                no3_flow_fraction * vno3mov * layer_depth_ratio,
                (1.0 - no3_flow_fraction) * vno3mov,
            ),
            (
                no3_flow_fraction * vnurmov * layer_depth_ratio,
                (1.0 - no3_flow_fraction) * vnurmov,
            ),
        )
    }

    fn transfer_mobile_nutrients_between_layers(
        legacy: &mut LegacyGlobalState,
        layer: usize,
        column: usize,
        old_water_content: f64,
        moved_water: f64,
        uplimit: f64,
        no3_flow_fraction: f64,
        layer_depth_ratio: f64,
    ) {
        if moved_water <= 0.0 {
            return;
        }
        let nitconc =
            nutrient_concentration(legacy.vol_no3_n_content[[layer, column]], old_water_content);
        let nurconc = nutrient_concentration(
            legacy.vol_urea_n_content[[layer, column]],
            old_water_content,
        );
        legacy.vol_no3_n_content[[layer, column]] = uplimit * nitconc;
        legacy.vol_urea_n_content[[layer, column]] = uplimit * nurconc;

        let ((no3_to_lower, no3_retain), (urea_to_lower, urea_retain)) = move_mobile_nutrients(
            moved_water,
            nitconc,
            nurconc,
            no3_flow_fraction,
            layer_depth_ratio,
        );
        legacy.vol_no3_n_content[[layer + 1, column]] += no3_to_lower;
        legacy.vol_no3_n_content[[layer, column]] += no3_retain;
        legacy.vol_urea_n_content[[layer + 1, column]] += urea_to_lower;
        legacy.vol_urea_n_content[[layer, column]] += urea_retain;
    }

    fn drain_bottom_cell(
        legacy: &mut LegacyGlobalState,
        layer: usize,
        column: usize,
        old_water_content: f64,
        layer_depth: f64,
        col_width: f64,
    ) -> f64 {
        if legacy.vol_water_content[[layer, column]] <= legacy.max_water_capacity[layer] {
            return 0.0;
        }
        let drainage = (legacy.vol_water_content[[layer, column]]
            - legacy.max_water_capacity[layer])
            * layer_depth
            * col_width;
        let nitconc =
            nutrient_concentration(legacy.vol_no3_n_content[[layer, column]], old_water_content);
        let nurconc = nutrient_concentration(
            legacy.vol_urea_n_content[[layer, column]],
            old_water_content,
        );
        let saved_nitrogen = (legacy.vol_no3_n_content[[layer, column]]
            + legacy.vol_urea_n_content[[layer, column]])
            * layer_depth
            * col_width;
        legacy.vol_water_content[[layer, column]] = legacy.max_water_capacity[layer];
        legacy.vol_no3_n_content[[layer, column]] = nitconc * legacy.max_water_capacity[layer];
        legacy.vol_urea_n_content[[layer, column]] = nurconc * legacy.max_water_capacity[layer];
        let remaining_nitrogen = (legacy.vol_no3_n_content[[layer, column]]
            + legacy.vol_urea_n_content[[layer, column]])
            * layer_depth
            * col_width;
        legacy.soil_nitrogen_loss += saved_nitrogen - remaining_nitrogen;
        drainage
    }

    let mut nlx = legacy.nl;
    if legacy.water_table_layer < nlx {
        nlx = legacy.water_table_layer;
    }
    if nlx <= 0 {
        return 0.0;
    }

    let num_cols = legacy.nk as usize;
    let mut oldvh2oc = vec![0.0; num_cols];
    oldvh2oc.copy_from_slice(
        legacy
            .vol_water_content
            .slice(s![nlx as usize - 1, 0..num_cols])
            .as_slice()
            .expect("water row should be contiguous"),
    );

    for l in 0..(nlx as usize).saturating_sub(1) {
        let upper_water_before = legacy
            .vol_water_content
            .slice(s![l, 0..num_cols])
            .to_owned();
        oldvh2oc.copy_from_slice(
            upper_water_before
                .as_slice()
                .expect("owned row should be contiguous"),
        );
        let avwl = upper_water_before
            .iter()
            .zip(legacy.wk.iter().take(num_cols))
            .map(|(&vol_water, &col_width)| vol_water * col_width / legacy.row_space)
            .sum::<f64>();

        let uplimit = legacy.max_water_capacity[l];
        let layer_depth_ratio = legacy.dl[l] / legacy.dl[l + 1];
        let no3_flow_fraction = legacy.no3_flow_fraction[l];
        if avwl > uplimit {
            let mut wmov = avwl - uplimit;
            wmov *= layer_depth_ratio;
            legacy
                .vol_water_content
                .slice_mut(s![l, 0..num_cols])
                .fill(uplimit);
            for k in 0..num_cols {
                let col_width = legacy.wk[k];
                legacy.vol_water_content[[l + 1, k]] +=
                    wmov * col_width * legacy.nk as f64 / legacy.row_space;

                let qvout = oldvh2oc[k] - uplimit;
                transfer_mobile_nutrients_between_layers(
                    legacy,
                    l,
                    k,
                    oldvh2oc[k],
                    qvout,
                    uplimit,
                    no3_flow_fraction,
                    layer_depth_ratio,
                );
            }
        } else {
            for k in 0..num_cols {
                if legacy.vol_water_content[[l, k]] > uplimit {
                    let wmov = legacy.vol_water_content[[l, k]] - uplimit;
                    legacy.vol_water_content[[l, k]] = uplimit;
                    legacy.vol_water_content[[l + 1, k]] += wmov * layer_depth_ratio;
                    transfer_mobile_nutrients_between_layers(
                        legacy,
                        l,
                        k,
                        oldvh2oc[k],
                        wmov,
                        uplimit,
                        no3_flow_fraction,
                        layer_depth_ratio,
                    );
                }
            }
        }
    }

    let mut drainage = 0.0;
    let bottom_layer = nlx as usize - 1;
    for k in 0..num_cols {
        drainage += drain_bottom_cell(
            legacy,
            bottom_layer,
            k,
            oldvh2oc[k],
            legacy.dl[bottom_layer],
            legacy.wk[k],
        );
    }

    drainage
}

/// Advances capillary and lateral water flow for the active soil grid.
pub fn capillary_flow() {
    let mut legacy = LegacyGlobalState::from_globals();
    let mut scratch = read_hydrology_scratch();
    let mut wk1 = [0.0; 40];
    if legacy.daynum <= legacy.day_start {
        scratch.capillary_numiter = 0;
        for layer in 0..legacy.nl as usize {
            wk1[layer] = 0.0;
        }
    }

    scratch.capillary_numiter += 1;
    let capillary_numiter = scratch.capillary_numiter;

    recompute_soil_psi_grid(&mut legacy);

    let mut nlx = legacy.nl;
    if legacy.water_table_layer < nlx {
        nlx = legacy.water_table_layer - 1;
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
    for layer in 0..legacy.nl as usize {
        dl1[layer] = legacy.dl[layer];
        thad1[layer] = legacy.thad[layer];
        thts1[layer] = legacy.thts[layer];
        ps1[layer] = legacy.pore_space[layer];
    }

    if nlx > 0 {
        let n_cols = legacy.nk as usize;
        for column in 0..n_cols {
            load_vertical_flux_profiles(
                &legacy,
                column,
                nlx as usize,
                &mut q1,
                &mut q01,
                &mut psi1,
                &mut nit,
                &mut nur,
            );

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
                legacy.noitr,
                &legacy.soil_horizon_num,
                &legacy.beta,
                &legacy.saturated_hyd_cond,
                &legacy.alpha,
                legacy.ratio_implicit,
                legacy.conmax,
                capillary_numiter,
            );
            nitrogen_flow(legacy.nl as usize, &q01, &q1, &dl1, &mut nit, &mut nur);

            store_vertical_flux_profiles(&mut legacy, column, nlx as usize, &q1, &psi1, &nit, &nur);
        }

        let mut pp1 = [0.0; 40];
        let mut qr1 = [0.0; 40];
        let mut qs1 = [0.0; 40];

        for layer in 0..nlx as usize {
            load_horizontal_flux_profiles(
                &legacy, layer, n_cols, &mut q1, &mut q01, &mut psi1, &mut nit, &mut nur, &mut qr1,
                &mut qs1, &mut pp1, &mut wk1,
            );

            water_flux(
                &mut q1,
                &mut psi1,
                &wk1,
                &qr1,
                &qs1,
                &pp1,
                n_cols,
                0,
                layer,
                legacy.noitr,
                &legacy.soil_horizon_num,
                &legacy.beta,
                &legacy.saturated_hyd_cond,
                &legacy.alpha,
                legacy.ratio_implicit,
                legacy.conmax,
                capillary_numiter,
            );
            nitrogen_flow(n_cols, &q01, &q1, &wk1, &mut nit, &mut nur);

            store_horizontal_flux_profiles(&mut legacy, layer, n_cols, &q1, &psi1, &nit, &nur);
        }
    }

    let water_drained_out = drain_with_legacy(&mut legacy);
    if water_drained_out > 0.0 {
        legacy.cum_water_drained += 10.0 * water_drained_out / legacy.row_space;
    }

    recompute_soil_psi_grid(&mut legacy);

    write_hydrology_scratch(scratch);
    legacy.write_to_globals();
}

fn recompute_soil_psi_grid(legacy: &mut LegacyGlobalState) {
    for layer in 0..legacy.nl as usize {
        let horizon = legacy.soil_horizon_num[layer] as usize;
        for column in 0..legacy.nk as usize {
            let water_content = legacy.vol_water_content[[layer, column]];
            legacy.soil_psi[[layer, column]] = psiq(
                water_content,
                legacy.thad[layer],
                legacy.thts[layer],
                legacy.alpha[horizon],
                legacy.beta[horizon],
            ) - PsiOsmotic(
                water_content,
                legacy.thts[layer],
                legacy.el_cond_sat_soil_today,
            );
        }
    }
}

fn load_vertical_flux_profiles(
    legacy: &LegacyGlobalState,
    column: usize,
    n_layers: usize,
    q1: &mut [f64; 40],
    q01: &mut [f64; 40],
    psi1: &mut [f64; 40],
    nit: &mut [f64; 40],
    nur: &mut [f64; 40],
) {
    for layer in 0..n_layers {
        q1[layer] = legacy.vol_water_content[[layer, column]];
        q01[layer] = legacy.vol_water_content[[layer, column]];
        psi1[layer] = legacy.soil_psi[[layer, column]]
            + PsiOsmotic(
                legacy.vol_water_content[[layer, column]],
                legacy.thts[layer],
                legacy.el_cond_sat_soil_today,
            );
        nit[layer] = legacy.vol_no3_n_content[[layer, column]];
        nur[layer] = legacy.vol_urea_n_content[[layer, column]];
    }
}

fn store_vertical_flux_profiles(
    legacy: &mut LegacyGlobalState,
    column: usize,
    n_layers: usize,
    q1: &[f64; 40],
    psi1: &[f64; 40],
    nit: &[f64; 40],
    nur: &[f64; 40],
) {
    for layer in 0..n_layers {
        legacy.vol_water_content[[layer, column]] = q1[layer];
        legacy.vol_no3_n_content[[layer, column]] = nit[layer];
        legacy.vol_urea_n_content[[layer, column]] = nur[layer];
        legacy.soil_psi[[layer, column]] = psi1[layer]
            - PsiOsmotic(
                legacy.vol_water_content[[layer, column]],
                legacy.thts[layer],
                legacy.el_cond_sat_soil_today,
            );
    }
}

#[allow(clippy::too_many_arguments)]
fn load_horizontal_flux_profiles(
    legacy: &LegacyGlobalState,
    layer: usize,
    n_cols: usize,
    q1: &mut [f64; 40],
    q01: &mut [f64; 40],
    psi1: &mut [f64; 40],
    nit: &mut [f64; 40],
    nur: &mut [f64; 40],
    qr1: &mut [f64; 40],
    qs1: &mut [f64; 40],
    pp1: &mut [f64; 40],
    wk1: &mut [f64; 40],
) {
    for column in 0..n_cols {
        q1[column] = legacy.vol_water_content[[layer, column]];
        q01[column] = legacy.vol_water_content[[layer, column]];
        psi1[column] = legacy.soil_psi[[layer, column]]
            + PsiOsmotic(
                legacy.vol_water_content[[layer, column]],
                legacy.thts[layer],
                legacy.el_cond_sat_soil_today,
            );
        qr1[column] = legacy.thad[layer];
        qs1[column] = legacy.thts[layer];
        pp1[column] = legacy.pore_space[layer];
        nit[column] = legacy.vol_no3_n_content[[layer, column]];
        nur[column] = legacy.vol_urea_n_content[[layer, column]];
        wk1[column] = legacy.wk[column];
    }
}

fn store_horizontal_flux_profiles(
    legacy: &mut LegacyGlobalState,
    layer: usize,
    n_cols: usize,
    q1: &[f64; 40],
    psi1: &[f64; 40],
    nit: &[f64; 40],
    nur: &[f64; 40],
) {
    for column in 0..n_cols {
        legacy.vol_water_content[[layer, column]] = q1[column];
        legacy.soil_psi[[layer, column]] = psi1[column]
            - PsiOsmotic(
                legacy.vol_water_content[[layer, column]],
                legacy.thts[layer],
                legacy.el_cond_sat_soil_today,
            );
        legacy.vol_no3_n_content[[layer, column]] = nit[column];
        legacy.vol_urea_n_content[[layer, column]] = nur[column];
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
/// Computes the amount of water requested by predicted irrigation rules.
pub fn ComputeIrrigation() {
    let target_stress = get_target_stress();
    if target_stress == -9999. {
        return;
    }

    let irrig_method = LegacyGlobalState::from_globals().irrig_method;
    if irrig_method == 2 {
        predict_drip_irrigation(target_stress);
    } else {
        predict_surface_irrigation(target_stress);
    }

    // If the amount of water to be applied (AppliedWater) is non zero update the date of last irrigation.
    let mut legacy = LegacyGlobalState::from_globals();
    if legacy.applied_water > 1e-5 {
        legacy.last_irrigation = legacy.daynum;
        legacy.write_to_globals();
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
/// Computes actual plant-root water uptake and updates soil potentials.
pub fn WaterUptake() {
    let mut legacy = LegacyGlobalState::from_globals();
    let num_layers = legacy.num_layers_with_roots as usize;
    let num_cols = legacy.nk as usize;
    let row_space = legacy.row_space;
    let light_intercept = legacy.light_intercept;
    let reference_transp = legacy.reference_transp;
    let average_soil_psi = legacy.average_soil_psi;
    let total_required_n = legacy.total_required_n;
    let per_plant_area = legacy.per_plant_area;
    let el_cond_sat_soil = legacy.el_cond_sat_soil_today;

    let mut dl_layer = vec![0.0f64; num_layers];
    let mut thad_layer = vec![0.0f64; num_layers];
    let mut thts_layer = vec![0.0f64; num_layers];
    let mut thetar_layer = vec![0.0f64; num_layers];
    let mut soil_horizon = vec![0usize; num_layers];
    let mut root_left = vec![0usize; num_layers];
    let mut root_right = vec![0usize; num_layers];
    for layer in 0..num_layers {
        dl_layer[layer] = legacy.dl[layer];
        thad_layer[layer] = legacy.thad[layer];
        thts_layer[layer] = legacy.thts[layer];
        thetar_layer[layer] = legacy.thetar[layer];
        soil_horizon[layer] = legacy.soil_horizon_num[layer] as usize;
        root_left[layer] = legacy.root_col_num_left[layer] as usize;
        root_right[layer] = legacy.root_col_num_right[layer] as usize;
    }

    let mut wk_col = vec![0.0f64; num_cols];
    for column in 0..num_cols {
        wk_col[column] = legacy.wk[column];
    }

    let alpha_horizon = legacy.alpha.clone();
    let beta_horizon = legacy.beta.clone();

    let mut root_uptake = Array2::<f64>::zeros((num_layers, num_cols));
    let mut vol_water = Array2::<f64>::zeros((num_layers, num_cols));
    let mut soil_psi = Array2::<f64>::zeros((num_layers, num_cols));
    let mut vol_no3 = Array2::<f64>::zeros((num_layers, num_cols));
    let mut vol_nh4 = Array2::<f64>::zeros((num_layers, num_cols));
    for layer in 0..num_layers {
        for column in 0..num_cols {
            root_uptake[[layer, column]] = legacy.root_wt_capbl_uptake[[layer, column]];
            vol_water[[layer, column]] = legacy.vol_water_content[[layer, column]];
            soil_psi[[layer, column]] = legacy.soil_psi[[layer, column]];
            vol_no3[[layer, column]] = legacy.vol_no3_n_content[[layer, column]];
            vol_nh4[[layer, column]] = legacy.vol_nh4_n_content[[layer, column]];
        }
    }

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

    for layer in 0..num_layers {
        for column in 0..num_cols {
            legacy.vol_water_content[[layer, column]] = vol_water[[layer, column]];
            legacy.soil_psi[[layer, column]] = soil_psi[[layer, column]];
            legacy.vol_no3_n_content[[layer, column]] = vol_no3[[layer, column]];
            legacy.vol_nh4_n_content[[layer, column]] = vol_nh4[[layer, column]];
        }
    }
    legacy.actual_transpiration = actual_transpiration;
    legacy.supply_no3_n = supply_no3;
    legacy.supply_nh4_n = supply_nh4;
    legacy.write_to_globals();
}
