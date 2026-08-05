//! Carbon allocation and organ-growth routines for a cotton plant.
//!
//! [`PlantGrowth`] coordinates potential and actual growth using atmospheric
//! conditions, management operations, root impedance, and legacy plant state.
//! The routines update shared model globals rather than returning a new plant
//! state; the daily engine synchronizes those writes before and after the plant
//! phase. Growth calculations are deterministic and perform no file I/O.

use crate::atmosphere::Atmosphere;
use crate::plant::root::{PotentialRootGrowth, RootGrowth};
use crate::utils::fmax;
use crate::LegacyGlobalState;
use crate::{
    nadj, pixda, pixdz, AbscisedLeafWeight, ActualBollGrowth, ActualBurrGrowth, ActualSquareGrowth,
    ActualStemGrowth, AdjAddHeightRate, AgeOfBoll, AgeOfPreFruNode, AgeOfSite, AirTemp,
    AvrgDailyTemp, BloomWeightLoss, BollWeight, BurrWeight, BurrWeightGreenBolls,
    BurrWeightOpenBolls, CarbonAllocatedForRootGrowth, CarbonStress, CottonWeightGreenBolls,
    CottonWeightOpenBolls, CumNetPhotosynth, DayEmerge, DayFirstDef, DayInc, DayTimeTemp, Daynum,
    DefoliantAppRate, DefoliationDate, DefoliationMethod, DensityFactor, ExtraCarbon,
    FruitFraction, FruitGrowthRatio, FruitingCode, GreenBollsLost, Kday, KdayAdjust, LeafAge,
    LeafArea, LeafAreaIndex, LeafAreaMainStem, LeafAreaNodes, LeafAreaPreFru, LeafWeightAreaRatio,
    LeafWeightMainStem, LeafWeightNodes, LeafWeightPreFru, LightIntercept, LwpMin, NStressFruiting,
    NStressRoots, NStressVeg, NetPhotosynthesis, NightTimeTemp, NodeLayer, NodeLayerPreFru,
    NumAdjustDays, NumFruitBranches, NumGreenBolls, NumNodes, NumOpenBolls, NumPreFruNodes,
    NumVegBranches, PerPlantArea, PercentDefoliation, PetioleWeightMainStem, PetioleWeightNodes,
    PetioleWeightPreFru, PlantHeight, PlantWeight, PlantWeightAtStart, PotGroAllBolls,
    PotGroAllBurrs, PotGroAllLeaves, PotGroAllPetioles, PotGroAllRoots, PotGroAllSquares,
    PotGroBolls, PotGroBurrs, PotGroLeafAreaMainStem, PotGroLeafAreaNodes, PotGroLeafAreaPreFru,
    PotGroLeafWeightMainStem, PotGroLeafWeightNodes, PotGroLeafWeightPreFru,
    PotGroPetioleWeightMainStem, PotGroPetioleWeightNodes, PotGroPetioleWeightPreFru,
    PotGroSquares, PotGroStem, ReserveC, RootWeightLoss, RowSpace, SquareWeight, StemWeight,
    TotalActualLeafGrowth, TotalActualPetioleGrowth, TotalLeafArea, TotalLeafWeight,
    TotalPetioleWeight, TotalRootWeight, TotalSquareWeight, TotalStemWeight, VarPar, WaterStress,
    WaterStressStem,
};
use cotton2k_core::RootImpedanceTables;

use super::Plant;
use crate::atmosphere::num_hours;
use crate::profile::AgronomyOperation;
use chrono::Duration;
use std::sync::{LazyLock, RwLock};

/// Ratio of actual leaf + petiole growth to their potential requirements.
/// Computed in `dry_matter_balance`, used by `actual_leaf_growth` (legacy name: `vratio`).
#[derive(Debug, Clone, Copy)]
struct GrowthScratch {
    v_ratio: f64,
    /// Defoliant intercepted by the canopy (kg/ha), accumulated across days.
    defkgh: f64,
    /// Cumulative defoliant amount used in the defoliation regression (kg/ha).
    tdfkgh: f64,
    /// Switch indicating whether a predicted defoliation date was set.
    idsw: i32,
}

impl Default for GrowthScratch {
    fn default() -> Self {
        Self {
            v_ratio: 1.0,
            defkgh: 0.0,
            tdfkgh: 0.0,
            idsw: 0,
        }
    }
}

static GROWTH_SCRATCH: LazyLock<RwLock<GrowthScratch>> =
    LazyLock::new(|| RwLock::new(GrowthScratch::default()));

fn read_growth_scratch() -> GrowthScratch {
    *GROWTH_SCRATCH
        .read()
        .expect("growth scratch state lock should not be poisoned")
}

fn with_growth_scratch_mut<R>(f: impl FnOnce(&mut GrowthScratch) -> R) -> R {
    let mut scratch = GROWTH_SCRATCH
        .write()
        .expect("growth scratch state lock should not be poisoned");
    f(&mut scratch)
}

pub(crate) fn reset_scratch_state() {
    with_growth_scratch_mut(|scratch| *scratch = GrowthScratch::default());
}

fn total_leaf_weight(legacy: &LegacyGlobalState) -> f64 {
    let mut result = 0.0;
    if legacy.first_square <= 0 {
        result += 0.2;
    }
    for i in 0..legacy.num_pre_fru_nodes.max(0) as usize {
        result += legacy.leaf_weight_pre_fru[i];
    }
    for k in 0..legacy.num_veg_branches.max(0) as usize {
        for l in 0..legacy.num_fruit_branches[k].max(0) as usize {
            result += legacy.leaf_weight_main_stem[[k, l]];
            for m in 0..legacy.num_nodes[[k, l]].max(0) as usize {
                result += legacy.leaf_weight_nodes[[k, l, m]];
            }
        }
    }
    result
}

fn total_leaf_area(legacy: &LegacyGlobalState) -> f64 {
    let mut result = 0.0;
    if legacy.first_square <= 0 {
        result += 0.12;
    }
    for i in 0..legacy.num_pre_fru_nodes.max(0) as usize {
        result += legacy.leaf_area_pre_fru[i];
    }
    for k in 0..legacy.num_veg_branches.max(0) as usize {
        for l in 0..legacy.num_fruit_branches[k].max(0) as usize {
            result += legacy.leaf_area_main_stem[[k, l]];
            for m in 0..legacy.num_nodes[[k, l]].max(0) as usize {
                result += legacy.leaf_area_nodes[[k, l, m]];
            }
        }
    }
    result
}

/// Integrates one day's potential and actual growth for a plant.
pub trait PlantGrowth {
    /// Simulate whole-plant growth for a day.
    fn grow(
        &mut self,
        atmosphere: Atmosphere,
        agronomy_ops: &[AgronomyOperation],
        root_tables: &RootImpedanceTables,
    );
    /// Applies a height increment and returns the resulting height change.
    fn plant_height_increment(&mut self, x: f64) -> f64;
}

impl PlantGrowth for Plant {
    /// This function simulates the potential and actual growth of cotton plants.
    /// It is called from the root engine's daily state transition and calls the following functions:
    /// actual_fruit_growth(), actual_leaf_growth(), ActualRootGrowth(),
    /// AddPlantHeight(), dry_matter_balance(), PotentialFruitGrowth(),
    /// potential_leaf_growth(), PotentialRootGrowth(), PotentialStemGrowth(), TotalLeafWeight().
    ///
    /// The following global variables are referenced here:
    /// ActualStemGrowth, DayInc, FirstSquare, FruitingCode, Kday, pixdz,
    /// PerPlantArea, RowSpace, WaterStressStem.
    ///
    /// The following global variables are set here:
    /// LeafAreaIndex, PlantHeight, PotGroAllRoots, PotGroStem, StemWeight,
    /// TotalLeafArea, TotalLeafWeight, TotalPetioleWeight, TotalStemWeight.
    fn grow(
        &mut self,
        atmosphere: Atmosphere,
        agronomy_ops: &[AgronomyOperation],
        root_tables: &RootImpedanceTables,
    ) {
        //     Call potential_leaf_growth() to compute potential growth rate of
        //     leaves.
        potential_leaf_growth();
        let mut legacy = LegacyGlobalState::from_globals();
        //     If it is after first square, call PotentialFruitGrowth() to compute
        //     potential
        //  growth rate of squares and bolls.
        if legacy.fruiting_code[[0, 0, 0]] > 0 {
            PotentialFruitGrowth(atmosphere.daylength);
            legacy = LegacyGlobalState::from_globals();
        }
        //     Active stem tissue (stemnew) is the difference between
        //     TotalStemWeight
        //  and the value of StemWeight(kkday).
        let voldstm = 32i32; // constant parameter (days for stem tissue to become "old")
        let mut kkday = legacy.kday - voldstm; // age of young stem tissue
        if kkday < 1 {
            kkday = 1;
        }
        let stemnew = legacy.total_stem_weight - legacy.stem_weight[kkday as usize]; // dry weight of active stem tissue.

        //     Call PotentialStemGrowth() to compute PotGroStem, potential growth
        //     rate of stems.
        //  The effect of temperature is introduced, by multiplying potential growth
        //  rate by DayInc. Stem growth is also affected by water stress
        //  (WaterStressStem) and possible PIX application (pixdz).   PotGroStem is
        //  limited by (maxstmgr * PerPlantArea) g per plant per day.
        legacy.pot_gro_stem =
            PotentialStemGrowth(stemnew) * legacy.day_inc * legacy.water_stress_stem * legacy.pixdz;
        let maxstmgr = 0.067; // maximum posible potential stem growth, g dm-2 day-1.
        if legacy.pot_gro_stem > maxstmgr * legacy.per_plant_area {
            legacy.pot_gro_stem = maxstmgr * legacy.per_plant_area;
        }
        legacy.write_to_globals();
        // Call PotentialRootGrowth() to compute potential growth rate of roots.
        // total potential growth rate of roots in g per slab. this
        // is computed in PotentialRootGrowth() and used in
        // ActualRootGrowth().
        let sumpdr = PotentialRootGrowth(root_tables);
        legacy = LegacyGlobalState::from_globals();
        // Total potential growth rate of roots is converted from g per slab (sumpdr) to g per plant (PotGroAllRoots).
        legacy.pot_gro_all_roots = sumpdr * 100. * legacy.per_plant_area / legacy.row_space;
        // Limit PotGroAllRoots to (maxrtgr*PerPlantArea) g per plant per day.
        let maxrtgr = 0.045; // maximum possible potential root growth, g dm-2 day-1.
        if legacy.pot_gro_all_roots > maxrtgr * legacy.per_plant_area {
            legacy.pot_gro_all_roots = maxrtgr * legacy.per_plant_area;
        }
        legacy.write_to_globals();
        // Call dry_matter_balance() to compute carbon balance, allocation of carbon to plant parts, and carbon stress. dry_matter_balance() also computes and returns the values of the following arguments:
        // cdleaf is carbohydrate requirement for leaf growth, g per plant per day. cdpet is carbohydrate requirement for petiole growth, g per plant per day. cdroot is carbohydrate requirement for root growth, g per plant per day. cdstem is carbohydrate requirement for stem growth, g per plant per day.
        let mut cdstem = 0.;
        let mut cdleaf = 0.;
        let mut cdpet = 0.;
        let mut cdroot = 0.;
        dry_matter_balance(&mut cdstem, &mut cdleaf, &mut cdpet, &mut cdroot);
        legacy = LegacyGlobalState::from_globals();
        // If it is after first square, call actual_fruit_growth() to compute actual growth rate of squares and bolls.
        if legacy.fruiting_code[[0, 0, 0]] > 0 {
            actual_fruit_growth();
            legacy = LegacyGlobalState::from_globals();
        }
        // Initialize TotalLeafWeight. It is assumed that cotyledons fall off at time of first square. Also initialize TotalLeafArea and TotalPetioleWeight.
        legacy.total_petiole_weight = 0.;
        legacy.write_to_globals();
        // Call actual_leaf_growth to compute actual growth rate of leaves and compute leaf area index.
        actual_leaf_growth();
        legacy = LegacyGlobalState::from_globals();
        legacy.leaf_area_index = total_leaf_area(&legacy) / legacy.per_plant_area;
        // Add ActualStemGrowth to TotalStemWeight, and define StemWeight(Kday) for this day.
        legacy.total_stem_weight += legacy.actual_stem_growth;
        legacy.stem_weight[legacy.kday as usize] = legacy.total_stem_weight;
        // Plant density affects growth in height of tall plants.
        let htdenf = 55.; // minimum plant height for plant density affecting growth in height.
                          // intermediate variable to compute denf2.
        let mut z1 = (legacy.plant_height - htdenf) / htdenf;
        if z1 < 0. {
            z1 = 0.;
        }
        if z1 > 1. {
            z1 = 1.;
        }
        // effect of plant density on plant growth in height.
        let denf2 = 1. + z1 * (legacy.density_factor - 1.);
        legacy.write_to_globals();
        // Call AddPlantHeight to compute PlantHeight.
        let height_inc = self.plant_height_increment(denf2);
        legacy = LegacyGlobalState::from_globals();
        legacy.plant_height += height_inc;
        legacy.write_to_globals();
        // Call ActualRootGrowth() to compute actual root growth.
        self.compute_actual_root_growth(sumpdr, agronomy_ops);
    }
    /// This function simulates the growth in height of the main stem of cotton plants.
    ///
    ///  It is called from PlantGrowth(). It returns the added plant height (cm).
    /// The following global variables are referenced here:
    ///  AdjAddHeightRate, AgeOfPreFruNode, AgeOfSite, CarbonStress, DayInc,
    ///  FruitingCode, Kday, KdayAdjust, nadj, NumAdjustDays, NumFruitBranches,
    ///  NumPreFruNodes, NStressVeg, pixdz, VarPar, WaterStressStem.
    ///  The argument used:
    ///   denf2 - effect of plant density on plant growth in height.
    fn plant_height_increment(&mut self, denf2: f64) -> f64 {
        //     The following constant parameters are used:
        const vhtpar: [f64; 7] = [1.0, 0.27, 0.60, 0.20, 0.10, 0.26, 0.32];
        let legacy = LegacyGlobalState::from_globals();
        let mut addz; // daily plant height growth increment, cm.
                      //     Calculate vertical growth of main stem before the square on the
                      //     second fruiting branch
                      //  has appeared. Added stem height (addz) is a function of the age of the
                      //  last prefruiting node.
        if legacy.fruiting_code[[0, 1, 0]] == 0 {
            addz = vhtpar[0]
                - vhtpar[1] * legacy.age_of_pre_fru_node[(legacy.num_pre_fru_nodes - 1) as usize];
            if addz > vhtpar[2] {
                addz = vhtpar[2];
            }
            if addz < 0. {
                addz = 0.;
            }
            //     It is assumed that the previous prefruiting node is also
            //  capable of growth, and its growth (dz2) is added to addz.
            if legacy.num_pre_fru_nodes > 1 {
                // plant height growth increment due to growth of the second node from the top.
                let dz2 = legacy.var_par[19]
                    - legacy.var_par[20]
                        * legacy.age_of_pre_fru_node[(legacy.num_pre_fru_nodes - 2) as usize];
                addz += if dz2 < 0. {
                    0.
                } else if dz2 > vhtpar[3] {
                    vhtpar[3]
                } else {
                    dz2
                };
            }
            //     The effect of water stress on stem height at this stage is
            //  less than at a later stage (as modified by vhtpar(4)).
            addz *= 1. - vhtpar[4] * (1. - legacy.water_stress_stem);
        } else {
            //     Calculate vertical growth of main stem after the second square
            //     has appeared.
            //  Added stem height (addz) is a function of the average  age (agetop)
            //  of the upper three main stem nodes.
            // node numbers of top three nodes.
            let l = (legacy.num_fruit_branches[0] - 1) as usize;
            let l1 = if l < 1 { 0 } else { l - 1 };
            let l2 = if l < 2 { 0 } else { l - 2 };
            // average physiological age of top three nodes.
            let agetop = (legacy.age_of_site[[0, l, 0]]
                + legacy.age_of_site[[0, l1, 0]]
                + legacy.age_of_site[[0, l2, 0]])
                / 3.;
            addz = legacy.var_par[21] + agetop * (legacy.var_par[22] + legacy.var_par[23] * agetop);
            if agetop > (-0.5 * legacy.var_par[22] / legacy.var_par[23]) {
                addz = legacy.var_par[24];
            }
            if addz < legacy.var_par[24] {
                addz = legacy.var_par[24];
            }
            if addz > legacy.var_par[25] {
                addz = legacy.var_par[25];
            }
            //     addz is affected by water, carbohydrate and nitrogen stresses.
            addz *= legacy.water_stress_stem;
            addz *= 1. - vhtpar[5] * (1. - legacy.carbon_stress);
            addz *= 1. - vhtpar[6] * (1. - legacy.n_stress_veg);
        }
        //     The effect of temperature is expressed by DayInc. there are also
        //     effects of
        //  pix, plant density, and of a variety-specific calibration parameter
        //  (VarPar(26)).
        addz *= legacy.var_par[26] * legacy.pixdz * legacy.day_inc * denf2;
        //    Apply adjustment to addz if plant map data have been read
        let kdadjustend = legacy.kday_adjust + legacy.num_adjust_days;
        if legacy.kday > legacy.kday_adjust && legacy.kday <= kdadjustend && legacy.nadj[1] {
            addz *= legacy.adj_add_height_rate;
        }
        addz
    }
}

/// Simulates potential leaf growth using a monomolecular leaf area curve.
///
/// Leaf area model (see `docs/plant-growth-variables.md`):
/// `L(t) = s_max * (1 - exp(-c * t^p))` and `r = dL/dt`.
/// This uses the temperature response `temperature_on_leaf_growth_rate`.
fn potential_leaf_growth() {
    const P: f64 = 1.6; // parameter of the leaf growth rate equation.
    const VPOTLF: [f64; 14] = [
        3.0, 0.95, 1.2, 13.5, -0.62143, 0.109365, 0.00137566, 0.025, 0.00005, 30., 0.02, 0.001,
        2.50, 0.18,
    ];
    let mut legacy = LegacyGlobalState::from_globals();
    let mut wstrlf =
        legacy.water_stress * (1. + VPOTLF[0] * (2. - legacy.water_stress)) - VPOTLF[0];
    if wstrlf < 0.05 {
        wstrlf = 0.05;
    }
    let wtfstrs = VPOTLF[1] + VPOTLF[2] * (1. - wstrlf);
    let mut tdday = legacy.avrg_daily_temp;
    if tdday < VPOTLF[3] {
        tdday = VPOTLF[3];
    }
    legacy.leaf_weight_area_ratio = wtfstrs / (VPOTLF[4] + tdday * (VPOTLF[5] - tdday * VPOTLF[6]));
    legacy.pot_gro_all_leaves = 0.;
    legacy.pot_gro_all_petioles = 0.;
    let mut c = 0.;
    let mut smax = 0.;
    for j in 0..legacy.num_pre_fru_nodes as usize {
        if legacy.leaf_area_pre_fru[j] <= 0. {
            legacy.pot_gro_leaf_area_pre_fru[j] = 0.;
            legacy.pot_gro_leaf_weight_pre_fru[j] = 0.;
            legacy.pot_gro_petiole_weight_pre_fru[j] = 0.;
        } else {
            let jp1 = (j + 1) as f64;
            smax = jp1 * (legacy.var_par[2] - legacy.var_par[3] * jp1);
            if smax < legacy.var_par[4] {
                smax = legacy.var_par[4];
            }
            c = VPOTLF[7] + VPOTLF[8] * jp1 * (jp1 - VPOTLF[9]);
            let age = legacy.age_of_pre_fru_node[j];
            let rate = smax * c * P * (-c * age.powf(P)).exp() * age.powf(P - 1.);
            if rate >= 1e-12 {
                legacy.pot_gro_leaf_area_pre_fru[j] = rate
                    * wstrlf
                    * legacy.pixda
                    * temperature_on_leaf_growth_rate(legacy.avrg_daily_temp);
                legacy.pot_gro_leaf_weight_pre_fru[j] =
                    legacy.pot_gro_leaf_area_pre_fru[j] * legacy.leaf_weight_area_ratio;
                legacy.pot_gro_petiole_weight_pre_fru[j] = legacy.pot_gro_leaf_area_pre_fru[j]
                    * legacy.leaf_weight_area_ratio
                    * VPOTLF[13];
                legacy.pot_gro_all_leaves += legacy.pot_gro_leaf_weight_pre_fru[j];
                legacy.pot_gro_all_petioles += legacy.pot_gro_petiole_weight_pre_fru[j];
            }
        }
    }
    let denfac = 1. - VPOTLF[12] * (1. - legacy.density_factor);
    for k in 0..legacy.num_veg_branches as usize {
        let nbrch = legacy.num_fruit_branches[k] as usize;
        for l in 0..nbrch {
            if legacy.leaf_area_main_stem[[k, l]] <= 0. {
                legacy.pot_gro_leaf_area_main_stem[[k, l]] = 0.;
                legacy.pot_gro_leaf_weight_main_stem[[k, l]] = 0.;
                legacy.pot_gro_petiole_weight_main_stem[[k, l]] = 0.;
            } else {
                let lp1 = (l + 1) as f64;
                smax = legacy.var_par[5] + legacy.var_par[6] * lp1 * (legacy.var_par[7] - lp1);
                smax *= denfac;
                if smax < legacy.var_par[4] {
                    smax = legacy.var_par[4];
                }
                c = VPOTLF[10] + lp1 * VPOTLF[11];
                let rate = if legacy.leaf_age[[k, l, 0]] > 70. {
                    0.
                } else {
                    let age = legacy.leaf_age[[k, l, 0]];
                    smax * c * P * (-c * age.powf(P)).exp() * age.powf(P - 1.)
                };
                if rate >= 1e-12 {
                    legacy.pot_gro_leaf_area_main_stem[[k, l]] = rate
                        * wstrlf
                        * legacy.pixda
                        * temperature_on_leaf_growth_rate(legacy.avrg_daily_temp);
                    legacy.pot_gro_leaf_weight_main_stem[[k, l]] =
                        legacy.pot_gro_leaf_area_main_stem[[k, l]] * legacy.leaf_weight_area_ratio;
                    legacy.pot_gro_petiole_weight_main_stem[[k, l]] = legacy
                        .pot_gro_leaf_area_main_stem[[k, l]]
                        * legacy.leaf_weight_area_ratio
                        * VPOTLF[13];
                    legacy.pot_gro_all_leaves += legacy.pot_gro_leaf_weight_main_stem[[k, l]];
                    legacy.pot_gro_all_petioles += legacy.pot_gro_petiole_weight_main_stem[[k, l]];
                }
            }
            let smaxx = smax;
            let cc = c;
            let nnid = legacy.num_nodes[[k, l]] as usize;
            for m in 0..nnid {
                if legacy.leaf_area_nodes[[k, l, m]] <= 0. {
                    legacy.pot_gro_leaf_area_nodes[[k, l, m]] = 0.;
                    legacy.pot_gro_leaf_weight_nodes[[k, l, m]] = 0.;
                    legacy.pot_gro_petiole_weight_nodes[[k, l, m]] = 0.;
                } else {
                    let mp1 = (m + 1) as f64;
                    smax = smaxx * (1. - legacy.var_par[8] * mp1);
                    c = cc * (1. - legacy.var_par[8] * mp1);
                    let rate = if legacy.leaf_age[[k, l, m]] > 70. {
                        0.
                    } else {
                        let age = legacy.leaf_age[[k, l, m]];
                        smax * c * P * (-c * age.powf(P)).exp() * age.powf(P - 1.)
                    };
                    if rate >= 1e-12 {
                        legacy.pot_gro_leaf_area_nodes[[k, l, m]] = rate
                            * wstrlf
                            * legacy.pixda
                            * temperature_on_leaf_growth_rate(legacy.avrg_daily_temp);
                        legacy.pot_gro_leaf_weight_nodes[[k, l, m]] = legacy
                            .pot_gro_leaf_area_nodes[[k, l, m]]
                            * legacy.leaf_weight_area_ratio;
                        legacy.pot_gro_petiole_weight_nodes[[k, l, m]] = legacy
                            .pot_gro_leaf_area_nodes[[k, l, m]]
                            * legacy.leaf_weight_area_ratio
                            * VPOTLF[13];
                        legacy.pot_gro_all_leaves += legacy.pot_gro_leaf_weight_nodes[[k, l, m]];
                        legacy.pot_gro_all_petioles +=
                            legacy.pot_gro_petiole_weight_nodes[[k, l, m]];
                    }
                }
            }
        }
    }
    legacy.write_to_globals();
}

/// Temperature response for leaf growth rate.
/// Returns a normalized multiplier in [0, 1] (see `docs/plant-growth-variables.md`).
fn temperature_on_leaf_growth_rate(t: f64) -> f64 {
    const PAR: [f64; 8] = [
        24.,
        -1.14277,
        0.0910026,
        0.00152344,
        -0.317136,
        0.0300712,
        0.000416356,
        0.2162044,
    ];
    let mut ra = if t > PAR[0] {
        PAR[1] + t * (PAR[2] - t * PAR[3])
    } else {
        PAR[4] + t * (PAR[5] - t * PAR[6])
    };
    if ra < 0. {
        ra = 0.;
    }
    ra / PAR[7]
}

/// Computes the cotton plant dry matter (carbon) balance and allocation.
///
/// Outputs include `CarbonStress`, organ allocations, `ExtraCarbon`,
/// `FruitGrowthRatio`, and `V_RATIO`. Formulas are documented in
/// `docs/plant-growth-variables.md`.
fn dry_matter_balance(cdstem: &mut f64, cdleaf: &mut f64, cdpet: &mut f64, cdroot: &mut f64) {
    let mut legacy = LegacyGlobalState::from_globals();
    const VCHBAL: [f64; 15] = [
        6.0, 2.5, 1.0, 5.0, 0.20, 0.80, 0.48, 0.40, 0.2072, 0.60651, 0.0065, 1.10, 4.0, 0.25, 4.0,
    ];
    let cdsqar =
        legacy.pot_gro_all_squares * (legacy.n_stress_fruiting + VCHBAL[0]) / (VCHBAL[0] + 1.);
    let cdboll = (legacy.pot_gro_all_bolls + legacy.pot_gro_all_burrs)
        * (legacy.n_stress_fruiting + VCHBAL[0])
        / (VCHBAL[0] + 1.);
    *cdleaf = legacy.pot_gro_all_leaves * (legacy.n_stress_veg + VCHBAL[1]) / (VCHBAL[1] + 1.);
    *cdstem = legacy.pot_gro_stem * (legacy.n_stress_veg + VCHBAL[2]) / (VCHBAL[2] + 1.);
    *cdroot = legacy.pot_gro_all_roots * (legacy.n_stress_roots + VCHBAL[3]) / (VCHBAL[3] + 1.);
    *cdpet = legacy.pot_gro_all_petioles * (legacy.n_stress_veg + VCHBAL[14]) / (VCHBAL[14] + 1.);
    let cdsum = *cdstem + *cdleaf + *cdpet + *cdroot + cdsqar + cdboll;
    let cpool = legacy.net_photosynthesis + legacy.reserve_c * VCHBAL[13];
    if cdsum <= 0. {
        legacy.carbon_stress = 1.;
        legacy.write_to_globals();
        return;
    }
    legacy.carbon_stress = cpool / cdsum;
    if legacy.carbon_stress > 1. {
        legacy.carbon_stress = 1.;
    }
    let mut pdboll = 0.;
    let mut pdsq = 0.;
    let mut xtrac1 = 0.;
    if legacy.carbon_stress >= 1. {
        legacy.total_actual_leaf_growth = *cdleaf;
        legacy.total_actual_petiole_growth = *cdpet;
        legacy.actual_stem_growth = *cdstem;
        legacy.carbon_allocated_for_root_growth = *cdroot;
        pdboll = cdboll;
        pdsq = cdsqar;
    } else {
        let mut cavail;
        if (cdboll + cdsqar) > 0. {
            let bsratio = cpool / (cdboll + cdsqar);
            let mut ffr = (VCHBAL[5] + VCHBAL[6] * (1. - legacy.water_stress)) * bsratio;
            if ffr < 0. {
                ffr = 0.;
            }
            if ffr > 1. {
                ffr = 1.;
            }
            if ffr > bsratio {
                ffr = bsratio;
            }
            pdboll = cdboll * ffr;
            pdsq = cdsqar * ffr;
            cavail = cpool - pdboll - pdsq;
        } else {
            cavail = cpool;
        }
        if (*cdleaf + *cdpet) > 0. {
            let mut flf = VCHBAL[7] * cavail / (*cdleaf + *cdpet);
            if flf < 0. {
                flf = 0.;
            }
            if flf > 1. {
                flf = 1.;
            }
            legacy.total_actual_leaf_growth = *cdleaf * flf;
            legacy.total_actual_petiole_growth = *cdpet * flf;
            cavail -= legacy.total_actual_leaf_growth + legacy.total_actual_petiole_growth;
        } else {
            legacy.total_actual_leaf_growth = 0.;
            legacy.total_actual_petiole_growth = 0.;
        }
        if *cdroot > 0. {
            let mut ratio = VCHBAL[8]
                + VCHBAL[9]
                    * (-VCHBAL[10]
                        * (legacy.total_stem_weight
                            + total_leaf_weight(&legacy)
                            + legacy.total_petiole_weight)
                        * legacy.per_plant_area)
                        .exp();
            ratio *= VCHBAL[11];
            let mut rtmax = ratio / (ratio + 1.);
            rtmax *= 1. + VCHBAL[12] * (1. - legacy.water_stress);
            if rtmax > 1. {
                rtmax = 1.;
            }
            let mut frt = rtmax * cavail / *cdroot;
            if frt < 0. {
                frt = 0.;
            }
            if frt > 1. {
                frt = 1.;
            }
            legacy.carbon_allocated_for_root_growth = fmax(*cdroot * frt, cavail - *cdstem);
            cavail -= legacy.carbon_allocated_for_root_growth;
        } else {
            legacy.carbon_allocated_for_root_growth = 0.;
        }
        if *cdstem > 0. {
            let mut fst = cavail / *cdstem;
            if fst < 0. {
                fst = 0.;
            }
            if fst > 1. {
                fst = 1.;
            }
            legacy.actual_stem_growth = *cdstem * fst;
        } else {
            legacy.actual_stem_growth = 0.;
        }
        if cavail > legacy.actual_stem_growth {
            xtrac1 = cavail - legacy.actual_stem_growth;
        }
    }
    if legacy.actual_stem_growth < 0. {
        legacy.actual_stem_growth = 0.;
    }
    if legacy.total_actual_leaf_growth < 0. {
        legacy.total_actual_leaf_growth = 0.;
    }
    if legacy.total_actual_petiole_growth < 0. {
        legacy.total_actual_petiole_growth = 0.;
    }
    if legacy.carbon_allocated_for_root_growth < 0. {
        legacy.carbon_allocated_for_root_growth = 0.;
    }
    if pdboll < 0. {
        pdboll = 0.;
    }
    if pdsq < 0. {
        pdsq = 0.;
    }
    legacy.reserve_c = legacy.reserve_c + legacy.net_photosynthesis
        - (legacy.actual_stem_growth
            + legacy.total_actual_leaf_growth
            + legacy.total_actual_petiole_growth
            + legacy.carbon_allocated_for_root_growth
            + pdboll
            + pdsq);
    let resmax = VCHBAL[4] * total_leaf_weight(&legacy);
    let xtrac2 = if legacy.reserve_c > resmax {
        let extra = legacy.reserve_c - resmax;
        legacy.reserve_c = resmax;
        extra
    } else {
        0.
    };
    legacy.extra_carbon = xtrac1 + xtrac2;
    if (legacy.pot_gro_all_squares + legacy.pot_gro_all_bolls + legacy.pot_gro_all_burrs) > 0. {
        legacy.fruit_growth_ratio = (pdsq + pdboll)
            / (legacy.pot_gro_all_squares + legacy.pot_gro_all_bolls + legacy.pot_gro_all_burrs);
    } else {
        legacy.fruit_growth_ratio = 1.;
    }
    with_growth_scratch_mut(|scratch| {
        if (legacy.pot_gro_all_leaves + legacy.pot_gro_all_petioles) > 0. {
            scratch.v_ratio = (legacy.total_actual_leaf_growth
                + legacy.total_actual_petiole_growth)
                / (legacy.pot_gro_all_leaves + legacy.pot_gro_all_petioles);
        } else {
            scratch.v_ratio = 1.;
        }
    });
    legacy.write_to_globals();
}

/// Simulates actual growth of squares and bolls.
/// Actual growth is proportional to potential growth via `FruitGrowthRatio`.
fn actual_fruit_growth() {
    let mut legacy = LegacyGlobalState::from_globals();
    legacy.total_square_weight = 0.;
    legacy.cotton_weight_green_bolls = 0.;
    legacy.burr_weight_green_bolls = 0.;
    legacy.actual_square_growth = 0.;
    legacy.actual_boll_growth = 0.;
    legacy.actual_burr_growth = 0.;
    for k in 0..legacy.num_veg_branches as usize {
        let nbrch = legacy.num_fruit_branches[k] as usize;
        for l in 0..nbrch {
            let nnid = legacy.num_nodes[[k, l]] as usize;
            for m in 0..nnid {
                if legacy.fruiting_code[[k, l, m]] == 1 {
                    let dwsq = legacy.pot_gro_squares[[k, l, m]] * legacy.fruit_growth_ratio;
                    legacy.square_weight[[k, l, m]] += dwsq;
                    legacy.actual_square_growth += dwsq;
                    legacy.total_square_weight += legacy.square_weight[[k, l, m]];
                }
                if legacy.fruiting_code[[k, l, m]] == 2 || legacy.fruiting_code[[k, l, m]] == 7 {
                    let dwboll = legacy.pot_gro_bolls[[k, l, m]] * legacy.fruit_growth_ratio;
                    legacy.boll_weight[[k, l, m]] += dwboll;
                    legacy.actual_boll_growth += dwboll;
                    legacy.cotton_weight_green_bolls += legacy.boll_weight[[k, l, m]];
                    let dwburr = legacy.pot_gro_burrs[[k, l, m]] * legacy.fruit_growth_ratio;
                    legacy.burr_weight[[k, l, m]] += dwburr;
                    legacy.actual_burr_growth += dwburr;
                    legacy.burr_weight_green_bolls += legacy.burr_weight[[k, l, m]];
                }
            }
        }
    }
    legacy.write_to_globals();
}

/// Simulates actual leaf and petiole growth.
/// Actual growth scales potential growth by `V_RATIO`.
fn actual_leaf_growth() {
    let mut legacy = LegacyGlobalState::from_globals();
    let v_ratio = read_growth_scratch().v_ratio;
    for j in 0..legacy.num_pre_fru_nodes as usize {
        legacy.leaf_weight_pre_fru[j] += legacy.pot_gro_leaf_weight_pre_fru[j] * v_ratio;
        legacy.petiole_weight_pre_fru[j] += legacy.pot_gro_petiole_weight_pre_fru[j] * v_ratio;
        legacy.total_petiole_weight += legacy.petiole_weight_pre_fru[j];
        legacy.leaf_area_pre_fru[j] += legacy.pot_gro_leaf_area_pre_fru[j] * v_ratio;
        legacy.leaf_area[legacy.node_layer_pre_fru[j] as usize] += legacy.leaf_area_pre_fru[j];
    }
    for k in 0..legacy.num_veg_branches as usize {
        let nbrch = legacy.num_fruit_branches[k] as usize;
        for l in 0..nbrch {
            legacy.leaf_weight_main_stem[[k, l]] +=
                legacy.pot_gro_leaf_weight_main_stem[[k, l]] * v_ratio;
            legacy.petiole_weight_main_stem[[k, l]] +=
                legacy.pot_gro_petiole_weight_main_stem[[k, l]] * v_ratio;
            legacy.total_petiole_weight += legacy.petiole_weight_main_stem[[k, l]];
            legacy.leaf_area_main_stem[[k, l]] +=
                legacy.pot_gro_leaf_area_main_stem[[k, l]] * v_ratio;
            legacy.leaf_area[legacy.node_layer[[k, l]] as usize] +=
                legacy.leaf_area_main_stem[[k, l]];
            let nnid = legacy.num_nodes[[k, l]] as usize;
            for m in 0..nnid {
                legacy.leaf_weight_nodes[[k, l, m]] +=
                    legacy.pot_gro_leaf_weight_nodes[[k, l, m]] * v_ratio;
                legacy.petiole_weight_nodes[[k, l, m]] +=
                    legacy.pot_gro_petiole_weight_nodes[[k, l, m]] * v_ratio;
                legacy.total_petiole_weight += legacy.petiole_weight_nodes[[k, l, m]];
                legacy.leaf_area_nodes[[k, l, m]] +=
                    legacy.pot_gro_leaf_area_nodes[[k, l, m]] * v_ratio;
                legacy.leaf_area[legacy.node_layer[[k, l]] as usize] +=
                    legacy.leaf_area_nodes[[k, l, m]];
            }
        }
    }
    legacy.write_to_globals();
}

/// Checks the dry matter balance for diagnostic purposes.
/// See `docs/plant-growth-variables.md` for the balance equations.
/// Checks and records the plant dry-matter balance for diagnostics.
pub fn check_dry_matter_balance() {
    let mut legacy = LegacyGlobalState::from_globals();
    let avail = legacy.plant_weight_at_start + legacy.cum_net_photosynth;
    legacy.plant_weight = legacy.total_root_weight
        + legacy.total_stem_weight
        + legacy.cotton_weight_green_bolls
        + legacy.burr_weight_green_bolls
        + total_leaf_weight(&legacy)
        + legacy.total_petiole_weight
        + legacy.total_square_weight
        + legacy.cotton_weight_open_bolls
        + legacy.burr_weight_open_bolls
        + legacy.reserve_c;
    let used = legacy.plant_weight
        + legacy.green_bolls_lost
        + legacy.abscised_leaf_weight
        + legacy.bloom_weight_loss
        + legacy.root_weight_loss;
    let _chobal = avail - used;
    legacy.write_to_globals();
}

/// Simulates the effects of defoliating chemicals.
/// Uses persistent state (`DEFKGH`, `TDFKGH`, `IDSW`) across days.
/// Applies the accumulated effects of defoliant treatments to plant state.
pub fn defoliate() {
    const P1: f64 = -50.0;
    const P2: f64 = 0.525;
    const P3: f64 = 7.06;
    const P4: f64 = 0.85;
    const P5: f64 = 2.48;
    const P6: f64 = 0.0374;
    const P7: f64 = 0.0020;
    let mut legacy = LegacyGlobalState::from_globals();
    let mut scratch = read_growth_scratch();
    if legacy.daynum <= legacy.day_emerge {
        scratch.tdfkgh = 0.;
        scratch.defkgh = 0.;
        scratch.idsw = 0;
    }
    for i in 0..5 {
        if legacy.num_open_bolls > 0. && legacy.defoliant_app_rate[i] <= -99.9 {
            let open_ratio = (100. * legacy.num_open_bolls
                / (legacy.num_open_bolls + legacy.num_green_bolls))
                as i32;
            if i == 0 && scratch.idsw == 0 {
                if (legacy.daynum >= legacy.defoliation_date[i] && legacy.defoliation_date[0] > 0)
                    || open_ratio > legacy.defoliation_method[i]
                {
                    scratch.idsw = 1;
                    legacy.defoliation_date[i] = legacy.daynum;
                    legacy.defoliant_app_rate[1] = -99.9;
                    if legacy.daynum < legacy.day_first_def || legacy.day_first_def <= 0 {
                        legacy.day_first_def = legacy.daynum;
                    }
                    legacy.defoliation_method[i] = 0;
                }
            }
            if i >= 1
                && legacy.daynum == (legacy.defoliation_date[i - 1] + 10)
                && legacy.leaf_area_index >= 0.2
            {
                legacy.defoliation_date[i] = legacy.daynum;
                if i < 4 {
                    legacy.defoliant_app_rate[i + 1] = -99.9;
                }
                legacy.defoliation_method[i] = 0;
            }
        }
        if legacy.daynum == legacy.defoliation_date[i] {
            if legacy.defoliant_app_rate[i] < -99. {
                scratch.tdfkgh = 2.5;
            } else {
                if legacy.defoliation_method[i] == 0 {
                    scratch.defkgh += legacy.defoliant_app_rate[i] * 0.95 * 1.12085 * 0.75;
                } else {
                    scratch.defkgh +=
                        legacy.defoliant_app_rate[i] * legacy.light_intercept * 1.12085 * 0.75;
                }
                scratch.tdfkgh += scratch.defkgh;
            }
        }
        if legacy.defoliation_date[i] > 0 && legacy.daynum > legacy.day_first_def {
            let dum = -legacy.lwp_min * 10.;
            legacy.percent_defoliation = P1
                + P2 * legacy.avrg_daily_temp
                + P3 * scratch.tdfkgh
                + P4 * (legacy.daynum - legacy.day_first_def) as f64
                + P5 * dum
                - P6 * dum * dum
                + P7 * legacy.avrg_daily_temp
                    * scratch.tdfkgh
                    * (legacy.daynum - legacy.day_first_def) as f64
                    * dum;
            if legacy.percent_defoliation < 0. {
                legacy.percent_defoliation = 0.;
            }
            let perdmax = 40.;
            if legacy.percent_defoliation > perdmax {
                legacy.percent_defoliation = perdmax;
            }
        }
    }
    with_growth_scratch_mut(|state| *state = scratch);
    legacy.write_to_globals();
}
/// Computes physiological age.
/// This function returns the daily 'physiological age' increment, based on hourly temperatures. It is called each day by SimulateThisDay().
/// The following global variable is used here:
///
/// AirTemp[] = array of hourly temperatures.
/// Computes the day's physiological-age increment from hourly air temperature.
pub fn PhysiologicalAge() -> f64 {
    // The following constant Parameters are used in this function:
    const p1: f64 = 12.; // threshold temperature, C
    const p2: f64 = 14.; // temperature, C, above p1, for one physiological day.
    const p3: f64 = 1.5; // maximum value of a physiological day.
                         // The threshold value is assumed to be 12 C (p1). One physiological day is equivalent to a day with an average temperature of 26 C, and therefore the heat units are divided by 14 (p2).
    let legacy = LegacyGlobalState::from_globals();
    // A linear relationship is assumed between temperature and heat unit accumulation in the range of 12 C (p1) to 33 C (p2*p3+p1).
    // The effect of temperatures higher than 33 C is assumed to be equivalent to that of 33 C.
    let mut dayfd = 0.; // the daily contribution to physiological age (return value).
    for ihr in 0..24 {
        let tfd = (legacy.air_temp[ihr] - p1) / p2; // the hourly contribution to physiological age.
        dayfd += if tfd < 0. {
            0.
        } else if tfd > p3 {
            p3
        } else {
            tfd
        };
    }
    dayfd / 24.
}
/// This function computes and returns the resistance of leaves of cotton
/// plants to transpiration. It is assumed to be a function of leaf age.
/// It is called from LeafWaterPotential().
///
/// The input argument (agel) is leaf age in physiological days.
/// Computes leaf transpiration resistance from physiological leaf age.
pub fn LeafResistance(agel: f64) -> f64 {
    // The following constant parameters are used:
    const afac: f64 = 160.; // factor used for computing leaf resistance.
    const agehi: f64 = 94.; // higher limit for leaf age.
    const agelo: f64 = 48.; // lower limit for leaf age.
    const rlmin: f64 = 0.5; // minimum leaf resistance.

    rlmin
        + if agel <= agelo {
            0.
        } else if agel >= agehi {
            (agehi - agelo).powi(2) / afac
        } else {
            (agel - agelo) * (2. * agehi - agelo - agel) / afac
        }
}
/// Simulates the potential growth of fruiting sites of cotton plants.
/// It is called from PlantGrowth(). It calls TemperatureOnFruitGrowthRate()
///
/// The following global variables are referenced here:
///       AgeOfBoll, AgeOfSite, FruitingCode, FruitFraction,
///       NumFruitBranches, NumNodes,  NumVegBranches, DayTimeTemp,
///       NightTimeTemp, VarPar, WaterStress.
/// The following global variables are set here:
///       PotGroAllBolls, PotGroAllBurrs, PotGroAllSquares, PotGroBolls,
///       PotGroBurrs, PotGroSquares.
/// References:
/// * https://doi.org/10.1016/0378-4290(79)90019-4
/// * https://agris.fao.org/agris-search/search.do?recordID=US9323856
fn PotentialFruitGrowth(daylength: Duration) {
    let mut legacy = LegacyGlobalState::from_globals();
    // The constant parameters used:
    const vpotfrt: [f64; 5] = [0.72, 0.30, 3.875, 0.125, 0.17];
    // Compute tfrt for the effect of temperature on boll and burr growth rates.
    // Function [TemperatureOnFruitGrowthRate()] is used (with parameters derived from GOSSYM), for day time and night time temperatures, weighted by day and night lengths.
    // the effect of temperature on rate of boll, burr or square growth.
    let tfrt = (num_hours(daylength) * TemperatureOnFruitGrowthRate(legacy.day_time_temp)
        + (24. - num_hours(daylength)) * TemperatureOnFruitGrowthRate(legacy.night_time_temp))
        / 24.;
    // Assign zero to sums of potential growth of squares, bolls and burrs.
    legacy.pot_gro_all_squares = 0.;
    legacy.pot_gro_all_bolls = 0.;
    legacy.pot_gro_all_burrs = 0.;
    // Assign values for the boll growth equation parameters.
    // These are cultivar - specific.
    // maximum boll growth period (physiological days).
    let agemax = legacy.var_par[9];
    // maximum rate of boll (seed and lint) growth,g per boll per physiological day.
    let rbmax = legacy.var_par[10];
    // maximum possible boll (seed and lint) weight, g per boll.
    let wbmax = legacy.var_par[11];
    for k in 0..legacy.num_veg_branches as usize {
        for l in 0..legacy.num_fruit_branches[k] as usize {
            for m in 0..legacy.num_nodes[[k, l]] as usize {
                // Calculate potential square growth for node (k,l,m).
                // Sum potential growth rates of squares as PotGroAllSquares.
                if legacy.fruiting_code[[k, l, m]] == 1 {
                    // ratesqr is the rate of square growth, g per square per day.
                    // The routine for this is derived from GOSSYM, and so are the parameters used.
                    let ratesqr = tfrt
                        * vpotfrt[3]
                        * (-vpotfrt[2] + vpotfrt[3] * legacy.age_of_site[[k, l, m]]).exp();
                    legacy.pot_gro_squares[[k, l, m]] = ratesqr * legacy.fruit_fraction[[k, l, m]];
                    legacy.pot_gro_all_squares += legacy.pot_gro_squares[[k, l, m]];
                }
                // Growth of seedcotton is simulated separately from the growth of burrs.
                //
                // The logistic function is used to simulate growth of seedcotton.
                //
                // The constants of this function for cultivar 'Acala-SJ2', are based on the data of Marani (1979);
                // they are derived from calibration for other cultivars agemax is the age of the boll (in physiological days after bloom) at the time when the boll growth rate is maximal.
                //
                // rbmax is the potential maximum rate of boll growth (g seeds plus lint dry weight per physiological day) at this age.
                //
                // wbmax is the maximum potential weight of seed plus lint (g dry weight per boll).
                //
                // The auxiliary variable pex is computed as:
                //
                //     pex = exp(-4 * rbmax * (t - agemax) / wbmax)
                //
                //  where t is the physiological age of the boll after bloom (=agebol).
                //
                // Boll weight (seed plus lint) at age T, according to the logistic function is:
                //
                //     wbol = wbmax / (1 + pex)
                //
                // and the potential boll growth rate at this age will be the
                // derivative of this function:
                //
                //     ratebol = 4 * rbmax * pex / (1. + pex)**2
                else if legacy.fruiting_code[[k, l, m]] == 2
                    || legacy.fruiting_code[[k, l, m]] == 7
                {
                    // pex is an intermediate variable to compute boll growth.
                    let pex =
                        (-4. * rbmax * (legacy.age_of_boll[[k, l, m]] - agemax) / wbmax).exp();
                    // ratebol is the rate of boll (seed and lint) growth, g per boll per day.
                    let ratebol = 4. * rbmax * pex / (1. + pex).powi(2) * tfrt;
                    // Potential growth rate of the burrs is assumed to be constant (vpotfrt[4] g dry weight per day) until the boll reaches its final volume.
                    // This occurs at the age of 22 physiological days in 'Acala-SJ2'.
                    // Both ratebol and ratebur are modified by temperature (tfrt) and ratebur is also affected by water stress (wfdb).

                    // rate of burr growth, g per boll per day.
                    let ratebur = if legacy.age_of_boll[[k, l, m]] >= 22. {
                        0.
                    } else {
                        // Compute wfdb for the effect of water stress on burr growth rate.
                        // wfdb is the effect of water stress on rate of burr growth.
                        let wfdb = vpotfrt[0] + vpotfrt[1] * legacy.water_stress;
                        vpotfrt[4]
                            * tfrt
                            * if wfdb < 0. {
                                0.
                            } else if wfdb > 1. {
                                1.
                            } else {
                                wfdb
                            }
                    };
                    // Potential boll (seeds and lint) growth rate (ratebol) and potential burr growth rate (ratebur) are multiplied by FruitFraction to compute PotGroBolls and PotGroBurrs for node (k,l,m).
                    legacy.pot_gro_bolls[[k, l, m]] = ratebol * legacy.fruit_fraction[[k, l, m]];
                    legacy.pot_gro_burrs[[k, l, m]] = ratebur * legacy.fruit_fraction[[k, l, m]];
                    // Sum potential growth rates of bolls and burrs as PotGroAllBolls and PotGroAllBurrs, respectively.
                    legacy.pot_gro_all_bolls += legacy.pot_gro_bolls[[k, l, m]];
                    legacy.pot_gro_all_burrs += legacy.pot_gro_burrs[[k, l, m]];
                }
                // If these are not green bolls, their potential growth is 0. End loop.
                else {
                    legacy.pot_gro_bolls[[k, l, m]] = 0.;
                    legacy.pot_gro_burrs[[k, l, m]] = 0.;
                }
            }
        }
    }
    legacy.write_to_globals();
}
/// Computes the effect of air temperature (t) on growth rate of bolls in cotton plants.
/// It is called from PotentialFruitGrowth().
///
/// Some values computed by this function:
///
/// | t (C)  |      tfr        |
/// |--------|-----------------|
/// | 12     | 0.              |
/// | 15     | 0.336           |
/// | 20     | 0.751           |
/// | 25     | 0.978           |
/// | 26     | 1.              |
/// | 28.5   | 1.024 (maximum) |
/// | 30     | 1.016           |
/// | 35     | 0.866           |
/// | 40     | 0.527           |
/// | 45     | 0.              |
fn TemperatureOnFruitGrowthRate(t: f64) -> f64 {
    const p1: f64 = -2.041;
    const p2: f64 = 0.215;
    const p3: f64 = 0.00377;
    let tfr = p1 + t * (p2 - p3 * t);

    if tfr < 0. {
        0.
    } else {
        tfr
    }
}

#[cfg(test)]
mod tests {
    use super::temperature_on_leaf_growth_rate;

    fn assert_close(actual: f64, expected: f64, eps: f64) {
        assert!(
            (actual - expected).abs() <= eps,
            "expected {expected}, got {actual}"
        );
    }

    #[test]
    fn temperature_on_leaf_growth_rate_reference_points() {
        let eps = 1e-3;
        assert_close(temperature_on_leaf_growth_rate(12.), 0.0, eps);
        assert_close(temperature_on_leaf_growth_rate(16.), 0.2655638, eps);
        assert_close(temperature_on_leaf_growth_rate(20.), 0.5446032, eps);
        assert_close(temperature_on_leaf_growth_rate(24.), 0.7620185, eps);
        assert_close(temperature_on_leaf_growth_rate(27.), 0.9422215, eps);
        assert_close(temperature_on_leaf_growth_rate(30.), 1.0000352, eps);
        assert_close(temperature_on_leaf_growth_rate(36.), 0.7351625, eps);
        assert_close(temperature_on_leaf_growth_rate(42.), 0.0, eps);
    }
}
/// Computes and returns the potential stem growth of cotton plants.
/// It is called from PlantGrowth().
///
/// The following argument is used:
/// * `stemnew` - dry weight of active stem tissue.
/// The following global variables are referenced here:
/// * [DensityFactor]
/// * [FruitingCode]
/// * [Kday]
/// * [VarPar]
fn PotentialStemGrowth(stemnew: f64) -> f64 {
    let legacy = LegacyGlobalState::from_globals();
    // There are two periods for computation of potential stem growth:
    // (1) Before the appearance of a square on the third fruiting branch.
    // Potential stem growth is a functon of plant age (Kday, days from emergence).
    if legacy.fruiting_code[[0, 2, 0]] == 0 {
        legacy.var_par[12] * (legacy.var_par[13] + legacy.var_par[14] * legacy.kday as f64)
    }
    // (2) After the appearance of a square on the third fruiting branch.
    // It is assumed that all stem tissue that is more than 32 days old is not active.
    // Potential stem growth is a function of active stem tissue weight (stemnew), and plant density (denfac).
    else {
        // effect of plant density on stem growth rate.
        let denfac = 1. - legacy.var_par[15] * (1. - legacy.density_factor);
        fmax(denfac, 0.2) * legacy.var_par[16] * (legacy.var_par[17] + legacy.var_par[18] * stemnew)
    }
}
