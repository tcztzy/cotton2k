//! Root impedance, initiation, growth, aging, and death routines.
//!
//! The root sub-model combines [`RootImpedanceTables`] with soil water, nitrate,
//! temperature, and mechanical constraints to update root cohorts in legacy
//! state. [`PotentialRootGrowth`] computes daily demand and [`RootGrowth`]
//! supplies the plant-growth integration; both mutate shared model state and do
//! not perform I/O.

//                   THE COTTON ROOT SUB-MODEL.
//     The following is a documentation of the root sub-model used in COTTON2K.
//     It is
//  derived from the principles of RHIZOS, as implemented in GOSSYM and in
//  GLYCIM, and from some principles of ROOTSIMU (Hoogenboom and Huck, 1986). It
//  is devised to be generally applicable, and may be used with root systems of
//  different crops by redefining the parameters, which are set here as
//  constants, and some of them are set in function InitializeRootData(). These
//  parameters are of course specific for the crop species, and perhaps also for
//  cultivars or cultivar groups.
//
//     This is a two-dimensional model and it may be used with soil cells of
//     different
//  sizes. The grid can be defined by the modeler. The maximum numbers of layers
//  and columns are given by the parameters maxl and maxk, respectively. These
//  are set to 40 and 20, in this version of COTTON2K. The grid is set at beginning.
//
//     The whole slab is being simulated. Thus, non-symmetrical processes (such
//     as
//  side-dressing of fertilizers or drip-irrigation) can be handled. The plant
//  is assumed to be situated at the center of the soil slab, or off-center for
//  skip-rows. Adjoining soil slabs are considered as mirror-images of each
//  other. Alternate-row drip systems (or any other agricultural input similarly
//  situated) are located at one edge of the slab.
//
//     The root mass in each cell is made up of NumRootAgeGroups classes, whose
//     number is to be
//  defined by the modeler. The maximum number of classes is 3 in this version
//  of COTTON2K.
//
//     The following functions account for root growth morphology:
//     TapRootGrowth() describes
//  growth of the taproot, and LateralRootGrowth() describes growth of the
//  lateral roots.
//
//     The calling sequence of the root submodel modules is as follows:
//     InitializeRootData() is called from ReadInput()
//  at the start of the simulation during profile initialization.
//     PotentialRootGrowth() and ActualRootGrowth() are called each day from
//     PlantGrowth(). PotentialRootGrowth() calls RootImpedance(),
//     soil_mechanic_resistance(), soil_air_on_root_growth(),
//  soil_nitrate_on_root_growth(), SoilTemOnRootGrowth(), soil_water_on_root_growth().
//     ActualRootGrowth() calls RedistRootNewGrowth(), TapRootGrowth(),
//     LateralRootGrowth(),
//  RootAging(), RootDeath(), root_summation().
use crate::model_state::for_each_cell;
use crate::LegacyGlobalState;
use crate::{
    cgind, dl, maxk, maxl, nk, nl, pixcon, rlat1, rlat2, wk, ActualRootGrowth, BulkDensity,
    CarbonAllocatedForRootGrowth, CumPlantNLoss, DailyRootLoss, DayEmerge, Daynum,
    DepthLastRootLayer, ExtraCarbon, LastTaprootLayer, LateralRootFlag, NumLayersWithRoots,
    NumRootAgeGroups, PerPlantArea, PixInPlants, PlantRowColumn, PlantRowLocation, PoreSpace,
    PotGroRoots, RootAge, RootColNumLeft, RootColNumRight, RootGroFactor, RootImpede, RootNConc,
    RootNitrogen, RootWeight, RootWeightLoss, RowSpace, SoilHorizonNum, SoilPsi, SoilTempDailyAvrg,
    TapRootLength, TotalRootWeight, VolNo3NContent, VolWaterContent,
};
use chrono::Datelike;
use ndarray::prelude::*;
use ndarray::{Array, Array2, Ix2};

use super::Plant;
use crate::profile::AgronomyOperation;
pub use cotton2k_core::RootImpedanceTables;

/// Calculates soil mechanical impedance to root growth, rtimpd(l,k), for all soil cells.
///
/// It is called from PotentialRootGrowth().
/// The impedance is a function of bulk density and water content in each soil soil cell.
/// No changes have been made in the original GOSSYM code.
///
/// Computes the soil mechanical impedance lookup table for each soil cell using the
/// impedance curves provided in [`RootImpedanceTables`].
fn RootImpedance(legacy: &mut LegacyGlobalState, tables: &RootImpedanceTables) {
    if tables.is_empty() {
        return;
    }
    let water_len = tables.water_content.len();
    let density_len = tables.bulk_density.len();
    if water_len == 0 || density_len == 0 {
        return;
    }
    let n_layers = legacy.nl as usize;
    let n_cols = legacy.nk as usize;
    let (vol_water_content, root_impede) = (&legacy.vol_water_content, &mut legacy.root_impede);
    for l in 0..n_layers {
        let j = legacy.soil_horizon_num[l] as usize;
        let bd_layer = legacy.bulk_density[j];
        let jj = tables
            .bulk_density
            .iter()
            .position(|&value| bd_layer <= value)
            .unwrap_or(density_len - 1);
        let j1 = jj.min(density_len - 1);
        let j0 = j1.saturating_sub(1);
        let bd_j0 = tables.bulk_density[j0];
        let bd_j1 = tables.bulk_density[j1];
        let water_row = vol_water_content.slice(s![l, 0..n_cols]);
        let mut impede_row = root_impede.slice_mut(s![l, 0..n_cols]);
        for (k, impede_cell) in impede_row.iter_mut().enumerate() {
            let vh2o = water_row[k] / bd_layer;
            let ik = tables
                .water_content
                .iter()
                .position(|&value| vh2o <= value)
                .unwrap_or(water_len - 1);
            let i1 = ik.min(water_len - 1);
            let i0 = i1.saturating_sub(1);
            if j1 == 0 {
                if i1 == 0 || vh2o <= tables.water_content[i1] {
                    *impede_cell = tables.impedance[j1][i1];
                } else {
                    let water_i0 = tables.water_content[i0];
                    let water_i1 = tables.water_content[i1];
                    let water_den = water_i1 - water_i0;
                    if water_den.abs() < f64::EPSILON {
                        *impede_cell = tables.impedance[j1][i1];
                    } else {
                        *impede_cell = tables.impedance[j1][i0]
                            - (tables.impedance[j1][i0] - tables.impedance[j1][i1])
                                * (vh2o - water_i0)
                                / water_den;
                    }
                }
            } else {
                let bd_den = bd_j0 - bd_j1;
                if i1 == 0 || vh2o <= tables.water_content[i1] {
                    if bd_den.abs() < f64::EPSILON {
                        *impede_cell = tables.impedance[j1][i1];
                    } else {
                        *impede_cell = tables.impedance[j0][i1]
                            - (tables.impedance[j0][i1] - tables.impedance[j1][i1])
                                * (bd_j0 - bd_layer)
                                / bd_den;
                    }
                } else {
                    let water_i0 = tables.water_content[i0];
                    let water_i1 = tables.water_content[i1];
                    let water_den = water_i1 - water_i0;
                    if water_den.abs() < f64::EPSILON || bd_den.abs() < f64::EPSILON {
                        *impede_cell = tables.impedance[j1][i1];
                    } else {
                        let temp1 = tables.impedance[j0][i1]
                            - (tables.impedance[j0][i1] - tables.impedance[j1][i1])
                                * (bd_j0 - bd_layer)
                                / bd_den;
                        let temp2 = tables.impedance[j0][i0]
                            - (tables.impedance[j0][i0] - tables.impedance[j1][i1])
                                * (bd_j0 - bd_layer)
                                / bd_den;
                        *impede_cell = temp2 + (temp1 - temp2) * (vh2o - water_i0) / water_den;
                    }
                }
            }
        }
    }
}

/// Calculates soil mechanical resistance of cell (l, k).
/// Uses the minimum `RootImpede` of adjacent cells and clamps the response.
/// See `docs/plant-growth-variables.md` for the exact equation.
/// It is called from [PotentialRootGrowth].
fn soil_mechanic_resistance(legacy: &LegacyGlobalState, l: i32, k: i32) -> f64 {
    const P1: f64 = 1.046;
    const P2: f64 = 0.034554;
    const P3: f64 = 0.5;
    let lp1 = if l == legacy.nl - 1 { l } else { l + 1 };
    let kp1 = if k == legacy.nk - 1 { k } else { k + 1 };
    let km1 = if k == 0 { 0 } else { k - 1 };
    let l = l as usize;
    let k = k as usize;
    let lp1 = lp1 as usize;
    let kp1 = kp1 as usize;
    let km1 = km1 as usize;
    let rtimpd0 = legacy.root_impede[[l, k]];
    let rtimpdkm1 = legacy.root_impede[[l, km1]];
    let rtimpdkp1 = legacy.root_impede[[l, kp1]];
    let rtimpdlp1 = legacy.root_impede[[lp1, k]];
    let mut rtimpdmin = rtimpd0;
    if rtimpdkm1 < rtimpdmin {
        rtimpdmin = rtimpdkm1;
    }
    if rtimpdkp1 < rtimpdmin {
        rtimpdmin = rtimpdkp1;
    }
    if rtimpdlp1 < rtimpdmin {
        rtimpdmin = rtimpdlp1;
    }
    let mut rtpct = P1 - P2 * rtimpdmin;
    if rtpct > 1. {
        rtpct = 1.;
    }
    if rtpct < P3 {
        rtpct = P3;
    }
    rtpct
}

/// Calculates reduction of potential root growth rate due to low oxygen content.
/// See `docs/plant-growth-variables.md` for the exact equation.
fn soil_air_on_root_growth(psislk: f64, pore_space: f64, vh2oclk: f64) -> f64 {
    const P1: f64 = 0.;
    const P2: f64 = 1.;
    const P3: f64 = 0.1;
    let mut rtrdo = if psislk > P1 { P2 } else { 1. };
    if vh2oclk >= pore_space {
        rtrdo = P3;
    }
    rtrdo
}

/// Calculates reduction of potential root growth rate due to low nitrate content.
/// See `docs/plant-growth-variables.md` for the exact equation.
fn soil_nitrate_on_root_growth(vno3clk: f64) -> f64 {
    const P1: f64 = 0.;
    const P2: f64 = 1.;
    if vno3clk < P1 {
        P2
    } else {
        1.
    }
}

/// Returns the effect of soil moisture on cotton root potential growth rate.
/// See `docs/plant-growth-variables.md` for the exact equation.
fn soil_water_on_root_growth(psislk: f64) -> f64 {
    const P1: f64 = 20.;
    const P2: f64 = 16.;
    let mut smf = ((P1 + psislk) / P2).powi(3);
    if smf < 0.02 {
        smf = 0.02;
    }
    if smf > 1. {
        smf = 1.;
    }
    smf
}

/// Initiates lateral root growth based on distance from the taproot tip.
/// See `docs/plant-growth-variables.md` for the exact logic.
fn initiate_lateral_roots() {
    let mut legacy = LegacyGlobalState::from_globals();
    initiate_lateral_roots_in_state(&mut legacy);
    legacy.write_to_globals();
}

fn initiate_lateral_roots_in_state(legacy: &mut LegacyGlobalState) {
    const DISTLR: f64 = 12.;
    if legacy.last_taproot_layer < 0 {
        return;
    }
    let mut sdl = legacy.tap_root_length - legacy.depth_last_root_layer;
    for l in (0..=legacy.last_taproot_layer as usize).rev() {
        sdl += legacy.dl[l];
        if sdl > DISTLR && legacy.lateral_root_flag[l] == 1 {
            legacy.lateral_root_flag[l] = 2;
        }
    }
}

/// Summarizes root data for output/plotting and updates [TotalRootWeight].
/// See `docs/plant-growth-variables.md` for the exact equation.
fn root_summation() {
    let mut legacy = LegacyGlobalState::from_globals();
    root_summation_in_state(&mut legacy);
    legacy.write_to_globals();
}

fn root_summation_in_state(legacy: &mut LegacyGlobalState) {
    let roots = legacy
        .root_weight
        .slice(s![0..legacy.nl as usize, 0..legacy.nk as usize, 0..3])
        .sum();
    legacy.total_root_weight = roots * 100. * legacy.per_plant_area / legacy.row_space;
}

fn root_growth_capable_weight(legacy: &LegacyGlobalState, layer: usize, column: usize) -> f64 {
    let num_age_groups = legacy.num_root_age_groups as usize;
    legacy
        .root_weight
        .slice(s![layer, column, 0..num_age_groups])
        .iter()
        .zip(legacy.cgind.iter().take(num_age_groups))
        .map(|(&root_weight, &growth_index)| root_weight * growth_index)
        .sum()
}

/// Calculates the potential root growth rate.
/// The return value is the sum of potential root growth rates for the whole slab (sumpdr).
/// It is called from PlantGrowth(). It calls: RootImpedance(),
/// soil_nitrate_on_root_growth(), soil_air_on_root_growth(), soil_mechanic_resistance(),
/// SoilTemOnRootGrowth() and soil_water_on_root_growth().
///
/// The following global variables are referenced here:
///
/// cgind, Date, Daynum, NumLayersWithRoots, NumRootAgeGroups, nk,
/// PerPlantArea, PoreSpace, RootAge, RootWeight. SoilPsi,
/// SoilTempDailyAvrg, VolNo3NContent, VolWaterContent.
///
/// The following global variables are set here:    PotGroRoots, RootGroFactor
pub fn PotentialRootGrowth(tables: &RootImpedanceTables) -> f64 {
    // The following constant parameter is used:
    // potential relative growth rate of the roots (g/g/day).
    const rgfac: f64 = 0.36;
    let mut legacy = LegacyGlobalState::from_globals();
    // Initialize to zero the PotGroRoots array.
    for_each_cell(
        legacy.num_layers_with_roots as usize,
        legacy.nk as usize,
        |l, k| legacy.pot_gro_roots[[l, k]] = 0.,
    );
    RootImpedance(&mut legacy, tables);
    let mut sumpdr = 0.; // sum of potential root growth rate for the whole slab
    for_each_cell(
        legacy.num_layers_with_roots as usize,
        legacy.nk as usize,
        |l, k| {
            // Check if this soil cell contains roots (if RootAge is greater than 0), and execute the following if this is true.
            if legacy.root_age[[l, k]] <= 0. {
                return;
            }
            // In each soil cell with roots, the root weight capable of growth rtwtcg is computed as the sum of RootWeight[l][k][i] * cgind[i] for all root classes.
            let rtwtcg = root_growth_capable_weight(&legacy, l, k);
            // Compute the temperature factor for root growth by calling function SoilTemOnRootGrowth() for this layer.
            // soil temperature, C, this day's average for this cell.
            let stday = legacy.soil_temp_daily_avrg[[l, k]] - 273.161;
            // effect of soil temperature on root growth.
            let temprg = SoilTemOnRootGrowth(stday);
            // Compute soil mechanical resistance for each soil cell by calling SoilMechanicResistance{},
            // the effect of soil aeration on root growth by calling SoilAirOnRootGrowth(),
            // and the effect of soil nitrate on root growth by calling SoilNitrateOnRootGrowth().
            // effect of soil mechanical resistance on root growth (returned from SoilMechanicResistance).
            let rtpct = soil_mechanic_resistance(&legacy, l as i32, k as i32);
            // effect of oxygen deficiency on root growth (returned from SoilAirOnRootGrowth).
            let rtrdo = soil_air_on_root_growth(
                legacy.soil_psi[[l, k]],
                legacy.pore_space[l],
                legacy.vol_water_content[[l, k]],
            );
            // effect of nitrate deficiency on root growth (returned from SoilNitrateOnRootGrowth).
            let rtrdn = soil_nitrate_on_root_growth(legacy.vol_no3_n_content[[l, k]]);
            // The root growth resistance factor RootGroFactor(l,k), which can take a value between 0 and 1, is computed as the minimum of these resistance factors.
            // It is further modified by multiplying it by the soil moisture function SoilWaterOnRootGrowth().
            // Potential root growth PotGroRoots(l,k) in each cell is computed as a product of rtwtcg, rgfac, the temperature function temprg, and RootGroFactor(l,k).
            // It is also multiplied by PerPlantArea / 19.6, for the effect of plant population density on root growth:
            // it is made comparable to a population of 5 plants per m in 38" rows.
            // The sum of the potential growth for the whole slab is computed as sumpdr.
            let mut minres = if rtpct < rtrdo { rtpct } else { rtrdo };
            if rtrdn < minres {
                minres = rtrdn;
            }
            let rtpsi = soil_water_on_root_growth(legacy.soil_psi[[l, k]]);
            legacy.root_gro_factor[[l, k]] = rtpsi * minres;
            legacy.pot_gro_roots[[l, k]] =
                rtwtcg * rgfac * temprg * legacy.root_gro_factor[[l, k]] * legacy.per_plant_area
                    / 19.6;
            sumpdr += legacy.pot_gro_roots[[l, k]];
        },
    );
    legacy.write_to_globals();
    sumpdr
}

/// Integrates actual root growth and root-management operations.
pub trait RootGrowth {
    /// Applies allocated root carbon and agronomy operations to root cohorts.
    fn compute_actual_root_growth(&mut self, sumpdr: f64, agronomy_ops: &[AgronomyOperation]);
}

impl RootGrowth for Plant {
    /// This function calculates the actual root growth rate.
    /// It is called from function PlantGrowth().
    /// It calls the following functions:  initiate_lateral_roots(), LateralRootGrowthLeft(), LateralRootGrowthRight(), RedistRootNewGrowth(), RootAging(), RootDeath(), root_summation(), TapRootGrowth().
    ///
    /// The following global variables are referenced here:
    /// CarbonAllocatedForRootGrowth, cgind, DayEmerge,
    /// Daynum, DepthLastRootLayer, dl, ExtraCarbon, LateralRootFlag,
    /// LastTaprootLayer, nk, nl, NumLayersWithRoots, NumRootAgeGroups,
    /// PerPlantArea, pixcon, PlantRowColumn, PotGroRoots, RootAge, RootNConc,
    /// RowSpace, TapRootLength, wk.
    ///
    /// The following global variables are set here:
    ///
    /// ActualRootGrowth, CumPlantNLoss, DailyRootLoss, PixInPlants, RootNitrogen, RootWeight, RootWeightLoss.
    ///
    /// The arguments are:
    /// * `sumpdr` - Sum of potential root growth rate for the whole slab.
    /// * `agronomy_ops` - Agronomy operations scheduled for the season. Only cultivation
    ///   entries are acted upon.
    fn compute_actual_root_growth(&mut self, sumpdr: f64, agronomy_ops: &[AgronomyOperation]) {
        // The following constant parameters are used:
        // The index for the relative partitioning of root mass produced by new growth to class i.
        const RootGrowthIndex: [f64; 3] = [1.0, 0.0, 0.0];
        // the threshold ratio of root mass capable of growth to soil cell volume (g/cm3); when this threshold is reached, a part of root growth in this cell may be extended to adjoining cells.
        const rtminc: f64 = 0.0000001;
        let mut legacy = LegacyGlobalState::from_globals();
        // Assign zero to pavail if this is the day of emergence.
        if legacy.daynum <= legacy.day_emerge {
            self.pavail = 0.;
        }
        let mut adwr1 = Array2::zeros([maxl as usize, maxk as usize]); // actual growth rate from roots existing in this soil cell.
                                                                       // Assign zero to the arrays of actual root growth rate.
        for_each_cell(legacy.nl as usize, legacy.nk as usize, |l, k| {
            legacy.actual_root_growth[[l, k]] = 0.
        });
        // The amount of carbon allocated for root growth is calculated from
        // CarbonAllocatedForRootGrowth, converted to g dry matter per slab, and
        // added to previously allocated carbon that has not been used for growth
        // (pavail). If there is no potential root growth, this will be stored in
        // pavail. Otherwise, zero is assigned to pavail.
        if sumpdr <= 0. {
            self.pavail += legacy.carbon_allocated_for_root_growth * 0.01 * legacy.row_space
                / legacy.per_plant_area;
            legacy.write_to_globals();
            return;
        }
        // actual growth factor (ratio of available C to potential growth).
        // The ratio of available C to potential root growth (actgf) is calculated.
        // pavail (if not zero) is used here, and zeroed after being used.
        let actgf = (self.pavail
            + legacy.carbon_allocated_for_root_growth * 0.01 * legacy.row_space
                / legacy.per_plant_area)
            / sumpdr;
        self.pavail = 0.;
        for_each_cell(
            legacy.num_layers_with_roots as usize,
            legacy.nk as usize,
            |l, k| {
                // adwr1(l,k), is proportional to the potential growth rate of roots in this cell.
                if legacy.root_age[[l, k]] > 0. {
                    adwr1[[l, k]] = legacy.pot_gro_roots[[l, k]] * actgf;
                }
            },
        );
        // If extra carbon is available, it is assumed to be added to the taproot.
        if legacy.extra_carbon > 0. {
            // available carbon for taproot growth, in g dry matter per slab.
            //  ExtraCarbon is converted to availt (g dry matter per slab).
            let availt = legacy.extra_carbon * 0.01 * legacy.row_space / legacy.per_plant_area;
            let mut sdl = legacy.tap_root_length - legacy.depth_last_root_layer;
            // distance from the tip of the taproot, cm.
            let mut tpwt = Array::<f64, Ix2>::zeros([maxl as usize, 2]); // proportionality factors for allocating added dry matter among taproot soil cells.
            let mut sumwt = 0.; // sum of the tpwt array.
                                // Extra Carbon (availt) is added to soil cells with roots in the columns immediately to
                                // the left and to the right of the location of the plant row.
            for l in (0..legacy.last_taproot_layer as usize + 1).rev() {
                // The weighting factors for allocating the carbon (tpwt) are proportional to the
                // volume of each soil cell and its distance (sdl) from the tip of the taproot.
                sdl += legacy.dl[l];
                tpwt[[l, 0]] = sdl * legacy.dl[l] * legacy.wk[legacy.plant_row_column as usize];
                tpwt[[l, 1]] = sdl * legacy.dl[l] * legacy.wk[legacy.plant_row_column as usize + 1];
                sumwt += tpwt[[l, 0]] + tpwt[[l, 1]];
            }
            // The proportional amount of mass is added to the mass of the last (inactive) root class in each soil cell.
            for l in 0..legacy.last_taproot_layer as usize + 1 {
                legacy.root_weight[[
                    l,
                    legacy.plant_row_column as usize,
                    legacy.num_root_age_groups as usize - 1,
                ]] += availt * tpwt[[l, 0]] / sumwt;
                legacy.root_weight[[
                    l,
                    legacy.plant_row_column as usize + 1,
                    legacy.num_root_age_groups as usize - 1,
                ]] += availt * tpwt[[l, 1]] / sumwt;
            }
        }
        // Check each cell if the ratio of root weight capable of growth to cell volume (rtconc)
        // exceeds the threshold rtminc, and call RedistRootNewGrowth() for this cell.
        // Otherwise, all new growth is contained in the same cell, and the actual growth in this
        // cell, ActualRootGrowth(l,k) will be equal to adwr1(l,k).
        for_each_cell(legacy.nl as usize, legacy.nk as usize, |l, k| {
            if legacy.root_age[[l, k]] <= 0. {
                return;
            }
            // ratio of root weight capable of growth to cell volume.
            let mut rtconc = root_growth_capable_weight(&legacy, l, k);
            rtconc /= legacy.dl[l] * legacy.wk[k];
            if rtconc > rtminc {
                RedistRootNewGrowth(&mut legacy, l, k, adwr1[[l, k]]);
            } else {
                legacy.actual_root_growth[[l, k]] += adwr1[[l, k]];
            }
        });
        // The new actual growth ActualRootGrowth(l,k) in each cell is partitioned among the root
        // classes in it in proportion to the parameters RootGrowthIndex(i), and the previous values
        // of RootWeight(k,l,i), and added to RootWeight(k,l,i).
        let mut sumind = 0.; // sum of the growth index for all classes in a cell.
        for i in 0..legacy.num_root_age_groups as usize {
            sumind += RootGrowthIndex[i];
        }
        for_each_cell(
            legacy.num_layers_with_roots as usize,
            legacy.nk as usize,
            |l, k| {
                if legacy.root_age[[l, k]] <= 0. {
                    return;
                }
                // sum of growth index multiplied by root weight, for all classes in a cell.
                let mut sumgr = 0.;
                for i in 0..legacy.num_root_age_groups as usize {
                    sumgr += RootGrowthIndex[i] * legacy.root_weight[[l, k, i]];
                }
                for i in 0..legacy.num_root_age_groups as usize {
                    if sumgr > 0. {
                        legacy.root_weight[[l, k, i]] += legacy.actual_root_growth[[l, k]]
                            * RootGrowthIndex[i]
                            * legacy.root_weight[[l, k, i]]
                            / sumgr;
                    } else {
                        legacy.root_weight[[l, k, i]] +=
                            legacy.actual_root_growth[[l, k]] * RootGrowthIndex[i] / sumind;
                    }
                }
            },
        );
        // Call function TapRootGrowth() for taproot elongation, if the taproot has not already
        // reached the bottom of the slab.
        if legacy.last_taproot_layer < legacy.nl - 1
            || legacy.tap_root_length < legacy.depth_last_root_layer
        {
            TapRootGrowth(&mut legacy);
        }
        // Call functions for growth of lateral roots
        initiate_lateral_roots_in_state(&mut legacy);
        for l in 0..legacy.last_taproot_layer as usize {
            if legacy.lateral_root_flag[l] == 2 {
                LateralRootGrowthLeft(&mut legacy, l);
                LateralRootGrowthRight(&mut legacy, l);
            }
        }
        // Initialize DailyRootLoss (weight of sloughed roots) for this day.
        legacy.daily_root_loss = 0.;
        for_each_cell(
            legacy.num_layers_with_roots as usize,
            legacy.nk as usize,
            |l, k| {
                // Check RootAge to determine if this soil cell contains roots,
                // and then compute root aging and root death by calling RootAging() and RootDeath()
                // for each soil cell with roots.
                if legacy.root_age[[l, k]] > 0. {
                    RootAgingInState(&mut legacy, l, k);
                    RootDeathInState(&mut legacy, l, k);
                }
            },
        );
        // Check if cultivation is executed in this day and apply the prescribed soil disturbance.
        for operation in agronomy_ops {
            if let AgronomyOperation::cultivation { date, depth } = operation {
                if date.ordinal() as i32 == legacy.daynum {
                    root_cultivation_in_state(&mut legacy, *depth);
                }
            }
        }
        // Convert DailyRootLoss to g per plant units and add it to RootWeightLoss.
        legacy.daily_root_loss =
            legacy.daily_root_loss * 100. * legacy.per_plant_area / legacy.row_space;
        legacy.root_weight_loss += legacy.daily_root_loss;
        // Adjust RootNitrogen (root N content) and PixInPlants (plant Pix content)
        // for loss by death of roots.
        legacy.root_nitrogen -= legacy.daily_root_loss * legacy.root_n_conc;
        legacy.cum_plant_n_loss += legacy.daily_root_loss * legacy.root_n_conc;
        legacy.pix_in_plants -= legacy.daily_root_loss * legacy.pixcon;
        // Call function RootSummation().
        root_summation_in_state(&mut legacy);
        legacy.write_to_globals();
    }
}

/// Simulates the removal of roots caused by a mechanical cultivation event.
///
/// The routine removes root mass from surface soil layers down to the requested
/// depth for all columns more than 15 cm away from the plant row, and it
/// accumulates the lost mass in `DailyRootLoss`.
fn root_cultivation(depth_cm: f64) {
    let mut legacy = LegacyGlobalState::from_globals();
    root_cultivation_in_state(&mut legacy, depth_cm);
    legacy.write_to_globals();
}

fn root_cultivation_in_state(legacy: &mut LegacyGlobalState, depth_cm: f64) {
    if depth_cm <= 0. {
        return;
    }
    let mut cultivated_layers: usize = 0;
    let mut accumulated_depth = 0.;
    for layer_index in 0..legacy.nl as usize {
        accumulated_depth += legacy.dl[layer_index];
        if accumulated_depth >= depth_cm {
            cultivated_layers = layer_index;
            break;
        }
    }
    if accumulated_depth < depth_cm {
        cultivated_layers = legacy.nl as usize;
    }
    if cultivated_layers == 0 {
        return;
    }
    let mut distance_from_edge = 0.;
    let num_age_groups = legacy.num_root_age_groups as usize;
    for column_index in 0..legacy.nk as usize {
        distance_from_edge += legacy.wk[column_index];
        if (distance_from_edge - legacy.plant_row_location).abs() >= 15. {
            let mut cultivated_root_slice = legacy.root_weight.slice_mut(s![
                0..cultivated_layers,
                column_index,
                0..num_age_groups
            ]);
            legacy.daily_root_loss += cultivated_root_slice.sum();
            cultivated_root_slice.fill(0.);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::LegacyGlobalState;
    use std::sync::{Mutex, OnceLock};

    fn global_root_state_lock() -> std::sync::MutexGuard<'static, ()> {
        static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
        LOCK.get_or_init(|| Mutex::new(())).lock().unwrap()
    }

    #[test]
    fn cultivation_ignores_non_positive_depth() {
        let _guard = global_root_state_lock();
        let original = LegacyGlobalState::from_globals();
        let mut legacy = original.clone();
        legacy.nl = 1;
        legacy.nk = 1;
        legacy.num_root_age_groups = 1;
        legacy.dl[0] = 10.;
        legacy.wk[0] = 10.;
        legacy.plant_row_location = 0.;
        legacy.root_weight[[0, 0, 0]] = 1.0;
        legacy.daily_root_loss = 0.;
        legacy.write_to_globals();

        root_cultivation(0.0);

        let after = LegacyGlobalState::from_globals();
        assert_eq!(after.daily_root_loss, 0.0);
        assert_eq!(after.root_weight[[0, 0, 0]], 1.0);

        original.write_to_globals();
    }

    #[test]
    fn cultivation_removes_roots_beyond_row() {
        let _guard = global_root_state_lock();
        let original = LegacyGlobalState::from_globals();
        let mut legacy = original.clone();
        legacy.nl = 2;
        legacy.nk = 2;
        legacy.num_root_age_groups = 1;
        legacy.dl[0] = 10.;
        legacy.dl[1] = 10.;
        legacy.wk[0] = 10.;
        legacy.wk[1] = 10.;
        legacy.plant_row_location = 0.;
        legacy.root_weight[[0, 0, 0]] = 1.5;
        legacy.root_weight[[0, 1, 0]] = 2.0;
        legacy.daily_root_loss = 0.;
        legacy.write_to_globals();

        root_cultivation(15.0);

        let after = LegacyGlobalState::from_globals();
        assert!(
            (after.daily_root_loss - 2.0).abs() < 1e-9,
            "expected loss of 2.0, got {}",
            after.daily_root_loss
        );
        assert_eq!(after.root_weight[[0, 0, 0]], 1.5);
        assert_eq!(after.root_weight[[0, 1, 0]], 0.0);

        original.write_to_globals();
    }
}

/// Computes the redistribution of new growth of roots into adjacent soil cells.
/// It is called from ActualRootGrowth().
///
/// Redistribution is affected by the factors rgfdn, rgfsd, rgfup.
/// and the values of RootGroFactor(l,k) in this soil cell and in the adjacent cells.
/// The values of ActualRootGrowth(l,k) for this and for the adjacent soil cells are computed.
/// The code of this module is based, with major changes, on the code of GOSSYM.
///
/// The following arguments are referenced:
/// * `addwt` - actual growth rate of roots in this soil cell.
/// * `k`, `l` - column and layer numbers.
///
/// The following global variables are referenced here:
/// * [dl]
/// * [nk]
/// * [nl]
/// * [PlantRowColumn]
/// * [RootGroFactor]
/// * [wk]
///
/// The following global variables are set here:
/// * [ActualRootGrowth]
/// * [DepthLastRootLayer]
/// * [LastTaprootLayer]
/// * [NumLayersWithRoots]
/// * [RootAge]
/// * [RootColNumLeft]
/// * [RootColNumRight]
/// * [TapRootLength]
fn RedistRootNewGrowth(legacy: &mut LegacyGlobalState, l: usize, k: usize, addwt: f64) {
    // The following constant parameters are used.
    // These are relative factors for root growth to adjoining cells, downwards, sideways, and upwards, respectively.
    // These factors are relative to the volume of the soil cell from which growth originates.
    const rgfdn: f64 = 900.;
    const rgfsd: f64 = 600.;
    const rgfup: f64 = 10.;
    // Set the number of layer above and below this layer, and the number of columns to the right and to the left of this column.
    // layer above and below layer l.
    let lp1 = if l == legacy.nl as usize - 1 {
        l
    } else {
        l + 1
    };
    let lm1 = if l == 0 { l } else { l - 1 };
    // column to the left and to the right of column k.
    let kp1 = if k == legacy.nk as usize - 1 {
        k
    } else {
        k + 1
    };
    let km1 = if k == 0 { k } else { k - 1 };
    // Compute proportionality factors (efac1, efacl, efacr, efacu, efacd) as the product of RootGroFactor and the geotropic factors in the respective soil cells.
    // Note that the geotropic factors are relative to the volume of the soil cell.
    // Compute the sum srwp of the proportionality factors.
    // product of RootGroFactor and geotropic factor for this cell.
    let efac1 = legacy.dl[l] * legacy.wk[k] * legacy.root_gro_factor[[l, k]];
    // as efac1 for the cell to the left of this cell.
    let efacl = rgfsd * legacy.root_gro_factor[[l, km1]];
    // as efac1 for the cell to the right of this cell.
    let efacr = rgfsd * legacy.root_gro_factor[[l, kp1]];
    // as efac1 for the cell above this cell.
    let efacu = rgfup * legacy.root_gro_factor[[lm1, k]];
    // as efac1 for the cell below this cell.
    let efacd = rgfdn * legacy.root_gro_factor[[lp1, k]];
    // sum of all efac values.
    let srwp = efac1 + efacl + efacr + efacu + efacd;
    // If srwp is very small, all the added weight will be in the
    // same soil soil cell, and execution of this function is ended.
    if srwp < 1e-10 {
        legacy.actual_root_growth[[l, k]] = addwt;
        return;
    }
    // Allocate the added dry matter to this and the adjoining soil cells in proportion to the EFAC factors.
    legacy.actual_root_growth[[l, k]] += addwt * efac1 / srwp;
    legacy.actual_root_growth[[l, km1]] += addwt * efacl / srwp;
    legacy.actual_root_growth[[l, kp1]] += addwt * efacr / srwp;
    legacy.actual_root_growth[[lm1, k]] += addwt * efacu / srwp;
    legacy.actual_root_growth[[lp1, k]] += addwt * efacd / srwp;
    // If roots are growing into new soil soil cells, initialize their RootAge to 0.01.
    if legacy.root_age[[l, km1]] == 0. {
        legacy.root_age[[l, km1]] = 0.01;
    }
    if legacy.root_age[[l, kp1]] == 0. {
        legacy.root_age[[l, kp1]] = 0.01;
    }
    if legacy.root_age[[lm1, k]] == 0. {
        legacy.root_age[[lm1, k]] = 0.01;
    }
    // If this new compartment is in a new layer with roots, also initialize its RootColNumLeft and RootColNumRight values.
    if legacy.root_age[[lp1, k]] == 0. && efacd > 0. {
        legacy.root_age[[lp1, k]] = 0.01;
        if legacy.root_col_num_left[lp1] == 0 || k < legacy.root_col_num_left[lp1] as usize {
            legacy.root_col_num_left[lp1] = k as i32;
        }
        if legacy.root_col_num_right[lp1] == 0 || k > legacy.root_col_num_right[lp1] as usize {
            legacy.root_col_num_right[lp1] = k as i32;
        }
    }
    // If this is in the location of the taproot, and the roots reach a new soil layer,
    // update the taproot parameters TapRootLength, DepthLastRootLayer, and LastTaprootLayer.
    if k == legacy.plant_row_column as usize || k == legacy.plant_row_column as usize + 1 {
        if lp1 > legacy.last_taproot_layer as usize && efacd > 0. {
            legacy.tap_root_length = legacy.depth_last_root_layer + 0.01;
            legacy.depth_last_root_layer += legacy.dl[lp1];
            legacy.last_taproot_layer = lp1 as i32;
        }
    }
    // Update NumLayersWithRoots, if necessary, and the values of RootColNumLeft and RootColNumRight for this layer.
    if legacy.num_layers_with_roots <= l as i32 && efacd > 0. {
        legacy.num_layers_with_roots = l as i32 + 1;
    }
    if km1 < legacy.root_col_num_left[l] as usize {
        legacy.root_col_num_left[l] = km1 as i32;
    }
    if kp1 > legacy.root_col_num_right[l] as usize {
        legacy.root_col_num_right[l] = kp1 as i32;
    }
}
/// Called from ActualRootGrowth(). It updates the variable celage(l,k) for the age of roots in each soil cell containing roots. When root age
/// reaches a threshold thtrn(i), a transformation of root tissue from class i to class i+1 occurs. The proportion transformed is trn(i).
///
/// It has been adapted from the code of GOSSYM, but the threshold
/// age for this process is based on the time from when the roots first
/// grew into each soil cell (whereas the time from emergence was used
/// in GOSSYM). Note: only 3 root age groups are assumed here.
///
/// The following global variable is referenced here: SoilTempDailyAvrg. The
/// following global variables are set here:        RootAge, RootWeight. The
/// arguments k, l - are column and layer numbers.
fn RootAging(l: usize, k: usize) {
    let mut legacy = LegacyGlobalState::from_globals();
    RootAgingInState(&mut legacy, l, k);
    legacy.write_to_globals();
}

fn RootAgingInState(legacy: &mut LegacyGlobalState, l: usize, k: usize) {
    //     The following constant parameters are used:
    const thtrn: [f64; 2] = [20.0, 40.0]; // the time threshold, from the initial
                                          // penetration of roots to a soil cell, after which some of the root
                                          // mass of class i may be transferred into the next class (i+1).
    const trn: [f64; 2] = [0.0060, 0.0050]; //  the daily proportion of this transfer.
                                            //
                                            // daily average soil temperature (c) of soil cell.
    let stday = legacy.soil_temp_daily_avrg[[l, k]] - 273.161;
    legacy.root_age[[l, k]] += SoilTemOnRootGrowth(stday);
    //
    for i in 0..2 {
        if legacy.root_age[[l, k]] > thtrn[i] {
            // root mass transferred from one class to the next.
            let xtr = trn[i] * legacy.root_weight[[l, k, i]];
            legacy.root_weight[[l, k, i + 1]] += xtr;
            legacy.root_weight[[l, k, i]] -= xtr;
        }
    }
}
/// Computes the death of root tissue in each soil cell containing roots.
///
/// When root age reaches a threshold thdth(i), a proportion dth(i) of the roots in class i dies.
/// The mass of dead roots is added to DailyRootLoss.
/// It has been adapted from GOSSYM, but the threshold age for this process is based on the time from when the roots first grew into each soil cell.
/// It is assumed that root death rate is greater in dry soil, for all root classes except class 1.
/// Root death rate is increased to the maximum value in soil saturated with water.
///
///    The following global variables are referenced here:
///      RootAge, PoreSpace, SoilPsi, VolWaterContent
///    The following global variables are set here:
///      RootWeight, DailyRootLoss
///    The arguments k, l - are column and layer numbers.
fn RootDeath(l: usize, k: usize) {
    let mut legacy = LegacyGlobalState::from_globals();
    RootDeathInState(&mut legacy, l, k);
    legacy.write_to_globals();
}

fn RootDeathInState(legacy: &mut LegacyGlobalState, l: usize, k: usize) {
    // The constant parameters are used:
    // a parameter in the equation for computing dthfac.
    const aa: f64 = 0.008;
    // the daily proportion of death of root tissue.
    const dth: [f64; 3] = [0.0001, 0.0002, 0.0001];
    // a parameter in the equation for computing dthfac.
    const dthmax: f64 = 0.10;
    // a parameter in the equation for computing dthfac.
    const psi0: f64 = -14.5;
    // the time threshold, from the initial penetration of roots to a soil cell, after which death of root tissue of class i may occur.
    const thdth: [f64; 3] = [30.0, 50.0, 100.0];
    for i in 0..3 {
        if legacy.root_age[[l, k]] > thdth[i] {
            // the computed proportion of roots dying in each
            // class.
            let mut dthfac = dth[i];
            if legacy.vol_water_content[[l, k]] >= legacy.pore_space[l] {
                dthfac = dthmax;
            } else {
                if i <= 1 && legacy.soil_psi[[l, k]] <= psi0 {
                    dthfac += aa * (psi0 - legacy.soil_psi[[l, k]]);
                }
                if dthfac > dthmax {
                    dthfac = dthmax;
                }
            }
            legacy.daily_root_loss += legacy.root_weight[[l, k, i]] * dthfac;
            legacy.root_weight[[l, k, i]] -= legacy.root_weight[[l, k, i]] * dthfac;
        }
    }
}

/// Computes the elongation of the taproot. It is called from ActualRootGrowth(). It calls SoilTemOnRootGrowth().
///
/// The following global variables are referenced here:
///      dl, nl, NumRootAgeGroups, PlantRowColumn, PoreSpace, RootGroFactor,
///      SoilTempDailyAvrg, VolWaterContent.
///    The following global variables are set here:
///      DepthLastRootLayer, LastTaprootLayer, NumLayersWithRoots, RootAge,
///      RootColNumLeft, RootColNumRight, RootWeight, TapRootLength.
fn TapRootGrowth(legacy: &mut LegacyGlobalState) {
    //     The following constant parameters are used:
    const p1: f64 = 0.10; // constant parameter.
    const rtapr: f64 = 4.; // potential growth rate of the taproot, cm/day.
                           //     It is assumed that taproot elongation takes place irrespective
                           //  of the supply of carbon to the roots. This elongation occurs in the
                           //  two columns of the slab where the plant is located.
                           //     Tap root elongation does not occur in water logged soil (water
                           //     table).
                           //
                           // the second column in which taproot growth occurs.
    let klocp1 = legacy.plant_row_column as usize + 1;
    if legacy.vol_water_content[[
        legacy.last_taproot_layer as usize,
        legacy.plant_row_column as usize,
    ]] >= legacy.pore_space[legacy.last_taproot_layer as usize]
        || legacy.vol_water_content[[legacy.last_taproot_layer as usize, klocp1]]
            >= legacy.pore_space[legacy.last_taproot_layer as usize]
    {
        return;
    }
    //     Average soil resistance (avres) is computed at the root tip.
    // avres = average value of RootGroFactor for the two soil cells at the tip
    // of the taproot.
    let avres = (legacy.root_gro_factor[[
        legacy.last_taproot_layer as usize,
        legacy.plant_row_column as usize,
    ]] + legacy.root_gro_factor[[legacy.last_taproot_layer as usize, klocp1]])
        / 2.;
    //     It is assumed that a linear empirical function of avres controls the
    //     rate of
    //  taproot elongation. The potential elongation rate of the taproot is also
    //  modified by soil temperature (SoilTemOnRootGrowth function), soil
    //  resistance, and soil moisture near the root tip.
    //     Actual growth is added to the taproot length TapRootLength.
    // daily average soil temperature (C) at root tip.
    let stday = 0.5
        * (legacy.soil_temp_daily_avrg[[
            legacy.last_taproot_layer as usize,
            legacy.plant_row_column as usize,
        ]] + legacy.soil_temp_daily_avrg[[legacy.last_taproot_layer as usize, klocp1]])
        - 273.161;
    // added taproot length, cm
    let addtaprt = rtapr * (1. - p1 + avres * p1) * SoilTemOnRootGrowth(stday);
    legacy.tap_root_length += addtaprt;
    //     DepthLastRootLayer, the depth (in cm) to the end of the last layer
    //     with
    //  roots, is used to check if the taproot reaches a new soil layer.
    //  When the new value of TapRootLength is greater than DepthLastRootLayer -
    //  it means that the roots penetrate to a new soil layer. In this case, and
    //  if this is not the last layer in the slab, the following is executed:
    //     LastTaprootLayer and DepthLastRootLayer are incremented. If this is a
    //     new layer with
    //  roots, NumLayersWithRoots is also redefined and two soil cells of the
    //  new layer are defined as containing roots (by initializing
    //  RootColNumLeft and RootColNumRight).
    if legacy.last_taproot_layer > legacy.nl - 2
        || legacy.tap_root_length <= legacy.depth_last_root_layer
    {
        return;
    }
    //     The following is executed when the taproot reaches a new soil layer.
    legacy.last_taproot_layer += 1;
    legacy.depth_last_root_layer += legacy.dl[legacy.last_taproot_layer as usize];
    if legacy.last_taproot_layer > legacy.num_layers_with_roots - 1 {
        legacy.num_layers_with_roots = legacy.last_taproot_layer + 1;
        if legacy.num_layers_with_roots > legacy.nl {
            legacy.num_layers_with_roots = legacy.nl;
        }
    }
    if legacy.root_col_num_left[legacy.last_taproot_layer as usize] == 0
        || legacy.root_col_num_left[legacy.last_taproot_layer as usize] > legacy.plant_row_column
    {
        legacy.root_col_num_left[legacy.last_taproot_layer as usize] = legacy.plant_row_column;
    }
    if legacy.root_col_num_right[legacy.last_taproot_layer as usize] == 0
        || legacy.root_col_num_right[legacy.last_taproot_layer as usize] < klocp1 as i32
    {
        legacy.root_col_num_right[legacy.last_taproot_layer as usize] = klocp1 as i32;
    }
    //     RootAge is initialized for these soil cells.
    legacy.root_age[[
        legacy.last_taproot_layer as usize,
        legacy.plant_row_column as usize,
    ]] = 0.01;
    legacy.root_age[[legacy.last_taproot_layer as usize, klocp1]] = 0.01;
    //     Some of the mass of class 1 roots is transferred downwards to
    //  the new cells. The transferred mass is proportional to 2 cm of
    //  layer width, but it is not more than half the existing mass in the
    //  last layer.
    for i in 0..legacy.num_root_age_groups as usize {
        // root mass transferred to the cell below when the
        // elongating taproot reaches a new soil layer.
        // first column
        let mut tran = legacy.root_weight[[
            legacy.last_taproot_layer as usize - 1,
            legacy.plant_row_column as usize,
            i,
        ]] * 2.
            / legacy.dl[legacy.last_taproot_layer as usize - 1];
        if tran
            > 0.5
                * legacy.root_weight[[
                    legacy.last_taproot_layer as usize - 1,
                    legacy.plant_row_column as usize,
                    i,
                ]]
        {
            tran = 0.5
                * legacy.root_weight[[
                    legacy.last_taproot_layer as usize - 1,
                    legacy.plant_row_column as usize,
                    i,
                ]];
        }
        legacy.root_weight[[
            legacy.last_taproot_layer as usize,
            legacy.plant_row_column as usize,
            i,
        ]] += tran;
        legacy.root_weight[[
            legacy.last_taproot_layer as usize - 1,
            legacy.plant_row_column as usize,
            i,
        ]] -= tran;
        // second column
        tran = legacy.root_weight[[legacy.last_taproot_layer as usize - 1, klocp1, i]] * 2.
            / legacy.dl[legacy.last_taproot_layer as usize - 1];
        if tran > 0.5 * legacy.root_weight[[legacy.last_taproot_layer as usize - 1, klocp1, i]] {
            tran = 0.5 * legacy.root_weight[[legacy.last_taproot_layer as usize - 1, klocp1, i]];
        }
        legacy.root_weight[[legacy.last_taproot_layer as usize, klocp1, i]] += tran;
        legacy.root_weight[[legacy.last_taproot_layer as usize - 1, klocp1, i]] -= tran;
    }
}
/// This function computes the elongation of the lateral roots
/// in a soil layer(l) to the left. It is called from ActualRootGrowth().
///
/// It calls function SoilTemOnRootGrowth().
///
/// The following global variables are referenced here:
/// * NumRootAgeGroups
/// * PlantRowColumn
/// * PoreSpace
/// * RootGroFactor
/// * SoilTempDailyAvrg
/// * VolWaterContent
/// * wk
///
/// The following global variables are set here:
/// * RootAge
/// * RootColNumLeft
/// * RootWeight
///
/// The argument used:
/// l - layer number in the soil slab.
fn LateralRootGrowthLeft(legacy: &mut LegacyGlobalState, l: usize) {
    //     The following constant parameters are used:
    const p1: f64 = 0.10; // constant parameter.
    const rlatr: f64 = 3.6; // potential growth rate of lateral roots, cm/day.
    const rtran: f64 = 0.2; // the ratio of root mass transferred to a new soil
                            // soil cell, when a lateral root grows into it.
                            //     On its initiation, lateral root length is assumed to be equal to the
                            //  width of a soil column soil cell at the location of the taproot.
    if legacy.rlat1[l] <= 0. {
        legacy.rlat1[l] = legacy.wk[legacy.plant_row_column as usize];
    }
    // daily average soil temperature (C) at root tip.
    let stday = legacy.soil_temp_daily_avrg[[l, legacy.plant_row_column as usize]] - 273.161;
    // the effect of soil temperature on root growth.
    let temprg = SoilTemOnRootGrowth(stday);
    //     Define the column with the tip of this lateral root (ktip)
    let mut ktip = 0usize; // column with the tips of the laterals to the left
    let mut sumwk = 0.; // summation of columns width
    for k in (0..legacy.plant_row_column as usize + 1).rev() {
        sumwk += legacy.wk[k];
        if sumwk >= legacy.rlat1[l] {
            ktip = k;
            break;
        }
    }
    //     Compute growth of the lateral root to the left. Potential
    //  growth rate (u) is modified by the soil temperature function,
    //  and the linearly modified effect of soil resistance (RootGroFactor).
    //     Lateral root elongation does not occur in water logged soil.
    if legacy.vol_water_content[[l, ktip]] < legacy.pore_space[l] {
        legacy.rlat1[l] += rlatr * temprg * (1. - p1 + legacy.root_gro_factor[[l, ktip]] * p1);
        //     If the lateral reaches a new soil soil cell: a proportion (tran)
        //     of
        //	mass of roots is transferred to the new soil cell.
        if legacy.rlat1[l] > sumwk && ktip > 0 {
            // column into which the tip of the lateral grows to
            // left.
            let newktip = ktip - 1;
            for i in 0..legacy.num_root_age_groups as usize {
                let tran = legacy.root_weight[[l, ktip, i]] * rtran;
                legacy.root_weight[[l, ktip, i]] -= tran;
                legacy.root_weight[[l, newktip, i]] += tran;
            }
            //     RootAge is initialized for this soil cell.
            //     RootColNumLeft of this layer is redefined.
            if legacy.root_age[[l, newktip]] == 0. {
                legacy.root_age[[l, newktip]] = 0.01;
            }
            if newktip < legacy.root_col_num_left[l] as usize {
                legacy.root_col_num_left[l] = newktip as i32;
            }
        }
    }
}
/// Computes the elongation of the lateral roots in a soil layer(l) to the right. It is called from ActualRootGrowth().
/// It calls function SoilTemOnRootGrowth().
///
/// The following global variables are referenced here:
/// nk, NumRootAgeGroups, PlantRowColumn, PoreSpace, RootGroFactor,
/// SoilTempDailyAvrg, VolWaterContent, wk.
///
/// The following global variables are set here:
/// RootAge, RootColNumRight, RootWeight.
///
/// The argument used:      l - layer number in the soil slab.
fn LateralRootGrowthRight(legacy: &mut LegacyGlobalState, l: usize) {
    //     The following constant parameters are used:
    const p1: f64 = 0.10; // constant parameter.
    const rlatr: f64 = 3.6; // potential growth rate of lateral roots, cm/day.
    const rtran: f64 = 0.2; // the ratio of root mass transferred to a new soil
                            // soil cell, when a lateral root grows into it.
                            //     On its initiation, lateral root length is assumed to be equal to the
                            //     width
                            //  of a soil column soil cell at the location of the taproot.
    let klocp1 = legacy.plant_row_column as usize + 1;
    if legacy.rlat2[l] <= 0. {
        legacy.rlat2[l] = legacy.wk[klocp1];
    }
    // daily average soil temperature (C) at root tip.
    let stday = legacy.soil_temp_daily_avrg[[l, klocp1]] - 273.161;
    // the effect of soil temperature on root growth.
    let temprg = SoilTemOnRootGrowth(stday);
    // define the column with the tip of this lateral root (ktip)
    let mut ktip = 0usize; // column with the tips of the laterals to the right
    let mut sumwk = 0.;
    for k in klocp1..legacy.nk as usize {
        sumwk += legacy.wk[k];
        if sumwk >= legacy.rlat2[l] {
            ktip = k;
            break;
        }
    }
    //     Compute growth of the lateral root to the right. Potential
    //  growth rate is modified by the soil temperature function,
    //  and the linearly modified effect of soil resistance (RootGroFactor).
    //     Lateral root elongation does not occur in water logged soil.
    if legacy.vol_water_content[[l, ktip]] < legacy.pore_space[l] {
        legacy.rlat2[l] += rlatr * temprg * (1. - p1 + legacy.root_gro_factor[[l, ktip]] * p1);
        //     If the lateral reaches a new soil soil cell: a proportion (tran)
        //     of
        //	mass of roots is transferred to the new soil cell.
        if legacy.rlat2[l] > sumwk && ktip < legacy.nk as usize - 1 {
            // column into which the tip of the lateral grows to left.
            let newktip = ktip + 1; // column into which the tip of the lateral grows to left.
            for i in 0..legacy.num_root_age_groups as usize {
                let tran = legacy.root_weight[[l, ktip, i]] * rtran;
                legacy.root_weight[[l, ktip, i]] -= tran;
                legacy.root_weight[[l, newktip, i]] += tran;
            }
            //     RootAge is initialized for this soil cell.
            //     RootColNumLeft of this layer is redefined.
            if legacy.root_age[[l, newktip]] == 0. {
                legacy.root_age[[l, newktip]] = 0.01;
            }
            if newktip > legacy.root_col_num_right[l] as usize {
                legacy.root_col_num_right[l] = newktip as i32;
            }
        }
    }
}
/// Called from [PotentialRootGrowth()], [TapRootGrowth()], [LateralRootGrowthLeft()] and [LateralRootGrowthRight()].
/// It computes the effects of soil temperature on the rate growth.
/// It is essentially based on the usage of GOSSYM, but relative values are computed here.
/// The computed value returned by this function is between 0 and 1.
///
/// It is assumed that maximum root growth occurs at or above 30 C, and no root growth occurs at or below 13.5 C.
/// A quadratic response to temperature between these limits is assumed.
///
/// The following argument is used:
/// * `t` - Soil temperature (C), daily average.
///
/// The parameters used are p1, p2, p3, with the following results:
/// t =      14    16    18    20    22    24    26    28    30
/// trf =  .053  .261  .443  .600  .731  .837  .917  .971  1.00
fn SoilTemOnRootGrowth(t: f64) -> f64 {
    const p1: f64 = -2.12;
    const p2: f64 = 0.2;
    const p3: f64 = -0.0032;
    if t >= 30. {
        1.
    } else {
        let trf = p1 + t * (p2 + p3 * t);
        if trf > 1. {
            1.
        } else if trf < 0. {
            0.
        } else {
            trf
        }
    }
}
