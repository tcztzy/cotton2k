use crate::atmosphere::Atmosphere;
use crate::plant::root::{PotentialRootGrowth, RootGrowth, RootImpedanceTables};
use crate::utils::fmax;
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

use super::Plant;
use crate::atmosphere::num_hours;
use crate::profile::AgronomyOperation;
use chrono::Duration;

/// Ratio of actual leaf + petiole growth to their potential requirements.
/// Computed in `dry_matter_balance`, used by `actual_leaf_growth` (legacy name: `vratio`).
static mut V_RATIO: f64 = 1.0;
/// Defoliant intercepted by the canopy (kg/ha), accumulated across days.
static mut DEFKGH: f64 = 0.0;
/// Cumulative defoliant amount used in the defoliation regression (kg/ha).
static mut TDFKGH: f64 = 0.0;
/// Switch indicating whether a predicted defoliation date was set.
static mut IDSW: i32 = 0;

pub trait PlantGrowth {
    /// Simulate whole-plant growth for a day.
    unsafe fn grow(
        &mut self,
        atmosphere: Atmosphere,
        agronomy_ops: &[AgronomyOperation],
        root_tables: &RootImpedanceTables,
    );
    unsafe fn plant_height_increment(&mut self, x: f64) -> f64;
}

impl PlantGrowth for Plant {
    /// This function simulates the potential and actual growth of cotton plants.
    /// It is called from [Profile::simulate_this_day()], and it calls the following functions:
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
    unsafe fn grow(
        &mut self,
        atmosphere: Atmosphere,
        agronomy_ops: &[AgronomyOperation],
        root_tables: &RootImpedanceTables,
    ) {
        //     Call potential_leaf_growth() to compute potential growth rate of
        //     leaves.
        potential_leaf_growth();
        //     If it is after first square, call PotentialFruitGrowth() to compute
        //     potential
        //  growth rate of squares and bolls.
        if FruitingCode[0][0][0] > 0 {
            PotentialFruitGrowth(atmosphere.daylength);
        }
        //     Active stem tissue (stemnew) is the difference between
        //     TotalStemWeight
        //  and the value of StemWeight(kkday).
        let voldstm = 32i32; // constant parameter (days for stem tissue to become "old")
        let mut kkday = Kday - voldstm; // age of young stem tissue
        if kkday < 1 {
            kkday = 1;
        }
        let stemnew = TotalStemWeight - StemWeight[kkday as usize]; // dry weight of active stem tissue.

        //     Call PotentialStemGrowth() to compute PotGroStem, potential growth
        //     rate of stems.
        //  The effect of temperature is introduced, by multiplying potential growth
        //  rate by DayInc. Stem growth is also affected by water stress
        //  (WaterStressStem) and possible PIX application (pixdz).   PotGroStem is
        //  limited by (maxstmgr * PerPlantArea) g per plant per day.
        PotGroStem = PotentialStemGrowth(stemnew) * DayInc * WaterStressStem * pixdz;
        let maxstmgr = 0.067; // maximum posible potential stem growth, g dm-2 day-1.
        if PotGroStem > maxstmgr * PerPlantArea {
            PotGroStem = maxstmgr * PerPlantArea;
        }
        // Call PotentialRootGrowth() to compute potential growth rate of roots.
        // total potential growth rate of roots in g per slab. this
        // is computed in PotentialRootGrowth() and used in
        // ActualRootGrowth().
        let sumpdr = PotentialRootGrowth(root_tables);
        // Total potential growth rate of roots is converted from g per slab (sumpdr) to g per plant (PotGroAllRoots).
        PotGroAllRoots = sumpdr * 100. * PerPlantArea / RowSpace;
        // Limit PotGroAllRoots to (maxrtgr*PerPlantArea) g per plant per day.
        let maxrtgr = 0.045; // maximum possible potential root growth, g dm-2 day-1.
        if PotGroAllRoots > maxrtgr * PerPlantArea {
            PotGroAllRoots = maxrtgr * PerPlantArea;
        }
        // Call dry_matter_balance() to compute carbon balance, allocation of carbon to plant parts, and carbon stress. dry_matter_balance() also computes and returns the values of the following arguments:
        // cdleaf is carbohydrate requirement for leaf growth, g per plant per day. cdpet is carbohydrate requirement for petiole growth, g per plant per day. cdroot is carbohydrate requirement for root growth, g per plant per day. cdstem is carbohydrate requirement for stem growth, g per plant per day.
        let mut cdstem = 0.;
        let mut cdleaf = 0.;
        let mut cdpet = 0.;
        let mut cdroot = 0.;
        dry_matter_balance(&mut cdstem, &mut cdleaf, &mut cdpet, &mut cdroot);
        // If it is after first square, call actual_fruit_growth() to compute actual growth rate of squares and bolls.
        if FruitingCode[0][0][0] > 0 {
            actual_fruit_growth();
        }
        // Initialize TotalLeafWeight. It is assumed that cotyledons fall off at time of first square. Also initialize TotalLeafArea and TotalPetioleWeight.
        TotalPetioleWeight = 0.;
        // Call actual_leaf_growth to compute actual growth rate of leaves and compute leaf area index.
        actual_leaf_growth();
        LeafAreaIndex = TotalLeafArea() / PerPlantArea;
        // Add ActualStemGrowth to TotalStemWeight, and define StemWeight(Kday) for this day.
        TotalStemWeight += ActualStemGrowth;
        StemWeight[Kday as usize] = TotalStemWeight;
        // Plant density affects growth in height of tall plants.
        let htdenf = 55.; // minimum plant height for plant density affecting growth in height.
                          // intermediate variable to compute denf2.
        let mut z1 = (PlantHeight - htdenf) / htdenf;
        if z1 < 0. {
            z1 = 0.;
        }
        if z1 > 1. {
            z1 = 1.;
        }
        // effect of plant density on plant growth in height.
        let denf2 = 1. + z1 * (DensityFactor - 1.);
        // Call AddPlantHeight to compute PlantHeight.
        PlantHeight += self.plant_height_increment(denf2);
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
    unsafe fn plant_height_increment(&mut self, denf2: f64) -> f64 {
        //     The following constant parameters are used:
        const vhtpar: [f64; 7] = [1.0, 0.27, 0.60, 0.20, 0.10, 0.26, 0.32];
        let mut addz; // daily plant height growth increment, cm.
                      //     Calculate vertical growth of main stem before the square on the
                      //     second fruiting branch
                      //  has appeared. Added stem height (addz) is a function of the age of the
                      //  last prefruiting node.
        if FruitingCode[0][1][0] == 0 {
            addz = vhtpar[0] - vhtpar[1] * AgeOfPreFruNode[(NumPreFruNodes - 1) as usize];
            if addz > vhtpar[2] {
                addz = vhtpar[2];
            }
            if addz < 0. {
                addz = 0.;
            }
            //     It is assumed that the previous prefruiting node is also
            //  capable of growth, and its growth (dz2) is added to addz.
            if NumPreFruNodes > 1 {
                // plant height growth increment due to growth of the second node from the top.
                let dz2 = VarPar[19] - VarPar[20] * AgeOfPreFruNode[(NumPreFruNodes - 2) as usize];
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
            addz *= 1. - vhtpar[4] * (1. - WaterStressStem);
        } else {
            //     Calculate vertical growth of main stem after the second square
            //     has appeared.
            //  Added stem height (addz) is a function of the average  age (agetop)
            //  of the upper three main stem nodes.
            // node numbers of top three nodes.
            let l = (NumFruitBranches[0] - 1) as usize;
            let l1 = if l < 1 { 0 } else { l - 1 };
            let l2 = if l < 2 { 0 } else { l - 2 };
            // average physiological age of top three nodes.
            let agetop = (AgeOfSite[0][l][0] + AgeOfSite[0][l1][0] + AgeOfSite[0][l2][0]) / 3.;
            addz = VarPar[21] + agetop * (VarPar[22] + VarPar[23] * agetop);
            if agetop > (-0.5 * VarPar[22] / VarPar[23]) {
                addz = VarPar[24];
            }
            if addz < VarPar[24] {
                addz = VarPar[24];
            }
            if addz > VarPar[25] {
                addz = VarPar[25];
            }
            //     addz is affected by water, carbohydrate and nitrogen stresses.
            addz *= WaterStressStem;
            addz *= 1. - vhtpar[5] * (1. - CarbonStress);
            addz *= 1. - vhtpar[6] * (1. - NStressVeg);
        }
        //     The effect of temperature is expressed by DayInc. there are also
        //     effects of
        //  pix, plant density, and of a variety-specific calibration parameter
        //  (VarPar(26)).
        addz *= VarPar[26] * pixdz * DayInc * denf2;
        //    Apply adjustment to addz if plant map data have been read
        let kdadjustend = KdayAdjust + NumAdjustDays;
        if Kday > KdayAdjust && Kday <= kdadjustend {
            if nadj[1] {
                addz *= AdjAddHeightRate;
            }
        }
        //
        return addz;
    }
}

/// Simulates potential leaf growth using a monomolecular leaf area curve.
///
/// Leaf area model (see `docs/plant-growth-variables.md`):
/// `L(t) = s_max * (1 - exp(-c * t^p))` and `r = dL/dt`.
/// This uses the temperature response `temperature_on_leaf_growth_rate`.
unsafe fn potential_leaf_growth() {
    const P: f64 = 1.6; // parameter of the leaf growth rate equation.
    const VPOTLF: [f64; 14] = [
        3.0, 0.95, 1.2, 13.5, -0.62143, 0.109365, 0.00137566, 0.025, 0.00005, 30., 0.02, 0.001,
        2.50, 0.18,
    ];
    let mut wstrlf = WaterStress * (1. + VPOTLF[0] * (2. - WaterStress)) - VPOTLF[0];
    if wstrlf < 0.05 {
        wstrlf = 0.05;
    }
    let wtfstrs = VPOTLF[1] + VPOTLF[2] * (1. - wstrlf);
    let mut tdday = AvrgDailyTemp;
    if tdday < VPOTLF[3] {
        tdday = VPOTLF[3];
    }
    LeafWeightAreaRatio = wtfstrs / (VPOTLF[4] + tdday * (VPOTLF[5] - tdday * VPOTLF[6]));
    PotGroAllLeaves = 0.;
    PotGroAllPetioles = 0.;
    let mut c = 0.;
    let mut smax = 0.;
    for j in 0..NumPreFruNodes as usize {
        if LeafAreaPreFru[j] <= 0. {
            PotGroLeafAreaPreFru[j] = 0.;
            PotGroLeafWeightPreFru[j] = 0.;
            PotGroPetioleWeightPreFru[j] = 0.;
        } else {
            let jp1 = (j + 1) as f64;
            smax = jp1 * (VarPar[2] - VarPar[3] * jp1);
            if smax < VarPar[4] {
                smax = VarPar[4];
            }
            c = VPOTLF[7] + VPOTLF[8] * jp1 * (jp1 - VPOTLF[9]);
            let age = AgeOfPreFruNode[j];
            let rate = smax * c * P * (-c * age.powf(P)).exp() * age.powf(P - 1.);
            if rate >= 1e-12 {
                PotGroLeafAreaPreFru[j] =
                    rate * wstrlf * pixda * temperature_on_leaf_growth_rate(AvrgDailyTemp);
                PotGroLeafWeightPreFru[j] = PotGroLeafAreaPreFru[j] * LeafWeightAreaRatio;
                PotGroPetioleWeightPreFru[j] =
                    PotGroLeafAreaPreFru[j] * LeafWeightAreaRatio * VPOTLF[13];
                PotGroAllLeaves += PotGroLeafWeightPreFru[j];
                PotGroAllPetioles += PotGroPetioleWeightPreFru[j];
            }
        }
    }
    let denfac = 1. - VPOTLF[12] * (1. - DensityFactor);
    for k in 0..NumVegBranches as usize {
        let nbrch = NumFruitBranches[k] as usize;
        for l in 0..nbrch {
            if LeafAreaMainStem[k][l] <= 0. {
                PotGroLeafAreaMainStem[k][l] = 0.;
                PotGroLeafWeightMainStem[k][l] = 0.;
                PotGroPetioleWeightMainStem[k][l] = 0.;
            } else {
                let lp1 = (l + 1) as f64;
                smax = VarPar[5] + VarPar[6] * lp1 * (VarPar[7] - lp1);
                smax *= denfac;
                if smax < VarPar[4] {
                    smax = VarPar[4];
                }
                c = VPOTLF[10] + lp1 * VPOTLF[11];
                let rate = if LeafAge[k][l][0] > 70. {
                    0.
                } else {
                    let age = LeafAge[k][l][0];
                    smax * c * P * (-c * age.powf(P)).exp() * age.powf(P - 1.)
                };
                if rate >= 1e-12 {
                    PotGroLeafAreaMainStem[k][l] =
                        rate * wstrlf * pixda * temperature_on_leaf_growth_rate(AvrgDailyTemp);
                    PotGroLeafWeightMainStem[k][l] =
                        PotGroLeafAreaMainStem[k][l] * LeafWeightAreaRatio;
                    PotGroPetioleWeightMainStem[k][l] =
                        PotGroLeafAreaMainStem[k][l] * LeafWeightAreaRatio * VPOTLF[13];
                    PotGroAllLeaves += PotGroLeafWeightMainStem[k][l];
                    PotGroAllPetioles += PotGroPetioleWeightMainStem[k][l];
                }
            }
            let smaxx = smax;
            let cc = c;
            let nnid = NumNodes[k][l] as usize;
            for m in 0..nnid {
                if LeafAreaNodes[k][l][m] <= 0. {
                    PotGroLeafAreaNodes[k][l][m] = 0.;
                    PotGroLeafWeightNodes[k][l][m] = 0.;
                    PotGroPetioleWeightNodes[k][l][m] = 0.;
                } else {
                    let mp1 = (m + 1) as f64;
                    smax = smaxx * (1. - VarPar[8] * mp1);
                    c = cc * (1. - VarPar[8] * mp1);
                    let rate = if LeafAge[k][l][m] > 70. {
                        0.
                    } else {
                        let age = LeafAge[k][l][m];
                        smax * c * P * (-c * age.powf(P)).exp() * age.powf(P - 1.)
                    };
                    if rate >= 1e-12 {
                        PotGroLeafAreaNodes[k][l][m] =
                            rate * wstrlf * pixda * temperature_on_leaf_growth_rate(AvrgDailyTemp);
                        PotGroLeafWeightNodes[k][l][m] =
                            PotGroLeafAreaNodes[k][l][m] * LeafWeightAreaRatio;
                        PotGroPetioleWeightNodes[k][l][m] =
                            PotGroLeafAreaNodes[k][l][m] * LeafWeightAreaRatio * VPOTLF[13];
                        PotGroAllLeaves += PotGroLeafWeightNodes[k][l][m];
                        PotGroAllPetioles += PotGroPetioleWeightNodes[k][l][m];
                    }
                }
            }
        }
    }
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
unsafe fn dry_matter_balance(
    cdstem: &mut f64,
    cdleaf: &mut f64,
    cdpet: &mut f64,
    cdroot: &mut f64,
) {
    const VCHBAL: [f64; 15] = [
        6.0, 2.5, 1.0, 5.0, 0.20, 0.80, 0.48, 0.40, 0.2072, 0.60651, 0.0065, 1.10, 4.0, 0.25, 4.0,
    ];
    let cdsqar = PotGroAllSquares * (NStressFruiting + VCHBAL[0]) / (VCHBAL[0] + 1.);
    let cdboll =
        (PotGroAllBolls + PotGroAllBurrs) * (NStressFruiting + VCHBAL[0]) / (VCHBAL[0] + 1.);
    *cdleaf = PotGroAllLeaves * (NStressVeg + VCHBAL[1]) / (VCHBAL[1] + 1.);
    *cdstem = PotGroStem * (NStressVeg + VCHBAL[2]) / (VCHBAL[2] + 1.);
    *cdroot = PotGroAllRoots * (NStressRoots + VCHBAL[3]) / (VCHBAL[3] + 1.);
    *cdpet = PotGroAllPetioles * (NStressVeg + VCHBAL[14]) / (VCHBAL[14] + 1.);
    let cdsum = *cdstem + *cdleaf + *cdpet + *cdroot + cdsqar + cdboll;
    let cpool = NetPhotosynthesis + ReserveC * VCHBAL[13];
    if cdsum <= 0. {
        CarbonStress = 1.;
        return;
    }
    CarbonStress = cpool / cdsum;
    if CarbonStress > 1. {
        CarbonStress = 1.;
    }
    let mut pdboll = 0.;
    let mut pdsq = 0.;
    let mut xtrac1 = 0.;
    if CarbonStress >= 1. {
        TotalActualLeafGrowth = *cdleaf;
        TotalActualPetioleGrowth = *cdpet;
        ActualStemGrowth = *cdstem;
        CarbonAllocatedForRootGrowth = *cdroot;
        pdboll = cdboll;
        pdsq = cdsqar;
    } else {
        let mut cavail;
        if (cdboll + cdsqar) > 0. {
            let bsratio = cpool / (cdboll + cdsqar);
            let mut ffr = (VCHBAL[5] + VCHBAL[6] * (1. - WaterStress)) * bsratio;
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
            TotalActualLeafGrowth = *cdleaf * flf;
            TotalActualPetioleGrowth = *cdpet * flf;
            cavail -= TotalActualLeafGrowth + TotalActualPetioleGrowth;
        } else {
            TotalActualLeafGrowth = 0.;
            TotalActualPetioleGrowth = 0.;
        }
        if *cdroot > 0. {
            let mut ratio = VCHBAL[8]
                + VCHBAL[9]
                    * (-VCHBAL[10]
                        * (TotalStemWeight + TotalLeafWeight() + TotalPetioleWeight)
                        * PerPlantArea)
                        .exp();
            ratio *= VCHBAL[11];
            let mut rtmax = ratio / (ratio + 1.);
            rtmax *= 1. + VCHBAL[12] * (1. - WaterStress);
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
            CarbonAllocatedForRootGrowth = fmax(*cdroot * frt, cavail - *cdstem);
            cavail -= CarbonAllocatedForRootGrowth;
        } else {
            CarbonAllocatedForRootGrowth = 0.;
        }
        if *cdstem > 0. {
            let mut fst = cavail / *cdstem;
            if fst < 0. {
                fst = 0.;
            }
            if fst > 1. {
                fst = 1.;
            }
            ActualStemGrowth = *cdstem * fst;
        } else {
            ActualStemGrowth = 0.;
        }
        if cavail > ActualStemGrowth {
            xtrac1 = cavail - ActualStemGrowth;
        }
    }
    if ActualStemGrowth < 0. {
        ActualStemGrowth = 0.;
    }
    if TotalActualLeafGrowth < 0. {
        TotalActualLeafGrowth = 0.;
    }
    if TotalActualPetioleGrowth < 0. {
        TotalActualPetioleGrowth = 0.;
    }
    if CarbonAllocatedForRootGrowth < 0. {
        CarbonAllocatedForRootGrowth = 0.;
    }
    if pdboll < 0. {
        pdboll = 0.;
    }
    if pdsq < 0. {
        pdsq = 0.;
    }
    ReserveC = ReserveC + NetPhotosynthesis
        - (ActualStemGrowth
            + TotalActualLeafGrowth
            + TotalActualPetioleGrowth
            + CarbonAllocatedForRootGrowth
            + pdboll
            + pdsq);
    let resmax = VCHBAL[4] * TotalLeafWeight();
    let xtrac2 = if ReserveC > resmax {
        let extra = ReserveC - resmax;
        ReserveC = resmax;
        extra
    } else {
        0.
    };
    ExtraCarbon = xtrac1 + xtrac2;
    if (PotGroAllSquares + PotGroAllBolls + PotGroAllBurrs) > 0. {
        FruitGrowthRatio = (pdsq + pdboll) / (PotGroAllSquares + PotGroAllBolls + PotGroAllBurrs);
    } else {
        FruitGrowthRatio = 1.;
    }
    if (PotGroAllLeaves + PotGroAllPetioles) > 0. {
        V_RATIO = (TotalActualLeafGrowth + TotalActualPetioleGrowth)
            / (PotGroAllLeaves + PotGroAllPetioles);
    } else {
        V_RATIO = 1.;
    }
}

/// Simulates actual growth of squares and bolls.
/// Actual growth is proportional to potential growth via `FruitGrowthRatio`.
unsafe fn actual_fruit_growth() {
    TotalSquareWeight = 0.;
    CottonWeightGreenBolls = 0.;
    BurrWeightGreenBolls = 0.;
    ActualSquareGrowth = 0.;
    ActualBollGrowth = 0.;
    ActualBurrGrowth = 0.;
    for k in 0..NumVegBranches as usize {
        let nbrch = NumFruitBranches[k] as usize;
        for l in 0..nbrch {
            let nnid = NumNodes[k][l] as usize;
            for m in 0..nnid {
                if FruitingCode[k][l][m] == 1 {
                    let dwsq = PotGroSquares[k][l][m] * FruitGrowthRatio;
                    SquareWeight[k][l][m] += dwsq;
                    ActualSquareGrowth += dwsq;
                    TotalSquareWeight += SquareWeight[k][l][m];
                }
                if FruitingCode[k][l][m] == 2 || FruitingCode[k][l][m] == 7 {
                    let dwboll = PotGroBolls[k][l][m] * FruitGrowthRatio;
                    BollWeight[k][l][m] += dwboll;
                    ActualBollGrowth += dwboll;
                    CottonWeightGreenBolls += BollWeight[k][l][m];
                    let dwburr = PotGroBurrs[k][l][m] * FruitGrowthRatio;
                    BurrWeight[k][l][m] += dwburr;
                    ActualBurrGrowth += dwburr;
                    BurrWeightGreenBolls += BurrWeight[k][l][m];
                }
            }
        }
    }
}

/// Simulates actual leaf and petiole growth.
/// Actual growth scales potential growth by `V_RATIO`.
unsafe fn actual_leaf_growth() {
    for j in 0..NumPreFruNodes as usize {
        LeafWeightPreFru[j] += PotGroLeafWeightPreFru[j] * V_RATIO;
        PetioleWeightPreFru[j] += PotGroPetioleWeightPreFru[j] * V_RATIO;
        TotalPetioleWeight += PetioleWeightPreFru[j];
        LeafAreaPreFru[j] += PotGroLeafAreaPreFru[j] * V_RATIO;
        LeafArea[NodeLayerPreFru[j] as usize] += LeafAreaPreFru[j];
    }
    for k in 0..NumVegBranches as usize {
        let nbrch = NumFruitBranches[k] as usize;
        for l in 0..nbrch {
            LeafWeightMainStem[k][l] += PotGroLeafWeightMainStem[k][l] * V_RATIO;
            PetioleWeightMainStem[k][l] += PotGroPetioleWeightMainStem[k][l] * V_RATIO;
            TotalPetioleWeight += PetioleWeightMainStem[k][l];
            LeafAreaMainStem[k][l] += PotGroLeafAreaMainStem[k][l] * V_RATIO;
            LeafArea[NodeLayer[k][l] as usize] += LeafAreaMainStem[k][l];
            let nnid = NumNodes[k][l] as usize;
            for m in 0..nnid {
                LeafWeightNodes[k][l][m] += PotGroLeafWeightNodes[k][l][m] * V_RATIO;
                PetioleWeightNodes[k][l][m] += PotGroPetioleWeightNodes[k][l][m] * V_RATIO;
                TotalPetioleWeight += PetioleWeightNodes[k][l][m];
                LeafAreaNodes[k][l][m] += PotGroLeafAreaNodes[k][l][m] * V_RATIO;
                LeafArea[NodeLayer[k][l] as usize] += LeafAreaNodes[k][l][m];
            }
        }
    }
}

/// Checks the dry matter balance for diagnostic purposes.
/// See `docs/plant-growth-variables.md` for the balance equations.
pub(crate) unsafe fn check_dry_matter_balance() {
    let avail = PlantWeightAtStart + CumNetPhotosynth;
    PlantWeight = TotalRootWeight
        + TotalStemWeight
        + CottonWeightGreenBolls
        + BurrWeightGreenBolls
        + TotalLeafWeight()
        + TotalPetioleWeight
        + TotalSquareWeight
        + CottonWeightOpenBolls
        + BurrWeightOpenBolls
        + ReserveC;
    let used = PlantWeight + GreenBollsLost + AbscisedLeafWeight + BloomWeightLoss + RootWeightLoss;
    let _chobal = avail - used;
}

/// Simulates the effects of defoliating chemicals.
/// Uses persistent state (`DEFKGH`, `TDFKGH`, `IDSW`) across days.
pub(crate) unsafe fn defoliate() {
    const P1: f64 = -50.0;
    const P2: f64 = 0.525;
    const P3: f64 = 7.06;
    const P4: f64 = 0.85;
    const P5: f64 = 2.48;
    const P6: f64 = 0.0374;
    const P7: f64 = 0.0020;
    if Daynum <= DayEmerge {
        TDFKGH = 0.;
        DEFKGH = 0.;
        IDSW = 0;
    }
    for i in 0..5 {
        if NumOpenBolls > 0. && DefoliantAppRate[i] <= -99.9 {
            let open_ratio = (100. * NumOpenBolls / (NumOpenBolls + NumGreenBolls)) as i32;
            if i == 0 && IDSW == 0 {
                if (Daynum >= DefoliationDate[i] && DefoliationDate[0] > 0)
                    || open_ratio > DefoliationMethod[i]
                {
                    IDSW = 1;
                    DefoliationDate[i] = Daynum;
                    DefoliantAppRate[1] = -99.9;
                    if Daynum < DayFirstDef || DayFirstDef <= 0 {
                        DayFirstDef = Daynum;
                    }
                    DefoliationMethod[i] = 0;
                }
            }
            if i >= 1 {
                if Daynum == (DefoliationDate[i - 1] + 10) && LeafAreaIndex >= 0.2 {
                    DefoliationDate[i] = Daynum;
                    if i < 4 {
                        DefoliantAppRate[i + 1] = -99.9;
                    }
                    DefoliationMethod[i] = 0;
                }
            }
        }
        if Daynum == DefoliationDate[i] {
            if DefoliantAppRate[i] < -99. {
                TDFKGH = 2.5;
            } else {
                if DefoliationMethod[i] == 0 {
                    DEFKGH += DefoliantAppRate[i] * 0.95 * 1.12085 * 0.75;
                } else {
                    DEFKGH += DefoliantAppRate[i] * LightIntercept * 1.12085 * 0.75;
                }
                TDFKGH += DEFKGH;
            }
        }
        if DefoliationDate[i] > 0 && Daynum > DayFirstDef {
            let dum = -LwpMin * 10.;
            PercentDefoliation = P1
                + P2 * AvrgDailyTemp
                + P3 * TDFKGH
                + P4 * (Daynum - DayFirstDef) as f64
                + P5 * dum
                - P6 * dum * dum
                + P7 * AvrgDailyTemp * TDFKGH * (Daynum - DayFirstDef) as f64 * dum;
            if PercentDefoliation < 0. {
                PercentDefoliation = 0.;
            }
            let perdmax = 40.;
            if PercentDefoliation > perdmax {
                PercentDefoliation = perdmax;
            }
        }
    }
}
/// Computes physiological age.
/// This function returns the daily 'physiological age' increment, based on hourly temperatures. It is called each day by SimulateThisDay().
/// The following global variable is used here:
///
/// AirTemp[] = array of hourly temperatures.
pub unsafe fn PhysiologicalAge() -> f64 {
    // The following constant Parameters are used in this function:
    const p1: f64 = 12.; // threshold temperature, C
    const p2: f64 = 14.; // temperature, C, above p1, for one physiological day.
    const p3: f64 = 1.5; // maximum value of a physiological day.
                         // The threshold value is assumed to be 12 C (p1). One physiological day is equivalent to a day with an average temperature of 26 C, and therefore the heat units are divided by 14 (p2).

    // A linear relationship is assumed between temperature and heat unit accumulation in the range of 12 C (p1) to 33 C (p2*p3+p1).
    // The effect of temperatures higher than 33 C is assumed to be equivalent to that of 33 C.
    let mut dayfd = 0.; // the daily contribution to physiological age (return value).
    for ihr in 0..24 {
        let tfd = (AirTemp[ihr] - p1) / p2; // the hourly contribution to physiological age.
        dayfd += if tfd < 0. {
            0.
        } else if tfd > p3 {
            p3
        } else {
            tfd
        };
    }
    return dayfd / 24.;
}
/// This function computes and returns the resistance of leaves of cotton
/// plants to transpiration. It is assumed to be a function of leaf age.
/// It is called from LeafWaterPotential().
///
/// The input argument (agel) is leaf age in physiological days.
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
unsafe fn PotentialFruitGrowth(daylength: Duration) {
    // The constant parameters used:
    const vpotfrt: [f64; 5] = [0.72, 0.30, 3.875, 0.125, 0.17];
    // Compute tfrt for the effect of temperature on boll and burr growth rates.
    // Function [TemperatureOnFruitGrowthRate()] is used (with parameters derived from GOSSYM), for day time and night time temperatures, weighted by day and night lengths.
    // the effect of temperature on rate of boll, burr or square growth.
    let tfrt = (num_hours(daylength) * TemperatureOnFruitGrowthRate(DayTimeTemp)
        + (24. - num_hours(daylength)) * TemperatureOnFruitGrowthRate(NightTimeTemp))
        / 24.;
    // Assign zero to sums of potential growth of squares, bolls and burrs.
    PotGroAllSquares = 0.;
    PotGroAllBolls = 0.;
    PotGroAllBurrs = 0.;
    // Assign values for the boll growth equation parameters.
    // These are cultivar - specific.
    // maximum boll growth period (physiological days).
    let agemax = VarPar[9];
    // maximum rate of boll (seed and lint) growth,g per boll per physiological day.
    let rbmax = VarPar[10];
    // maximum possible boll (seed and lint) weight, g per boll.
    let wbmax = VarPar[11];
    for k in 0..NumVegBranches as usize {
        for l in 0..NumFruitBranches[k] as usize {
            for m in 0..NumNodes[k][l] as usize {
                // Calculate potential square growth for node (k,l,m).
                // Sum potential growth rates of squares as PotGroAllSquares.
                if FruitingCode[k][l][m] == 1 {
                    // ratesqr is the rate of square growth, g per square per day.
                    // The routine for this is derived from GOSSYM, and so are the parameters used.
                    let ratesqr =
                        tfrt * vpotfrt[3] * (-vpotfrt[2] + vpotfrt[3] * AgeOfSite[k][l][m]).exp();
                    PotGroSquares[k][l][m] = ratesqr * FruitFraction[k][l][m];
                    PotGroAllSquares += PotGroSquares[k][l][m];
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
                else if FruitingCode[k][l][m] == 2 || FruitingCode[k][l][m] == 7 {
                    // pex is an intermediate variable to compute boll growth.
                    let pex = (-4. * rbmax * (AgeOfBoll[k][l][m] - agemax) / wbmax).exp();
                    // ratebol is the rate of boll (seed and lint) growth, g per boll per day.
                    let ratebol = 4. * rbmax * pex / (1. + pex).powi(2) * tfrt;
                    // Potential growth rate of the burrs is assumed to be constant (vpotfrt[4] g dry weight per day) until the boll reaches its final volume.
                    // This occurs at the age of 22 physiological days in 'Acala-SJ2'.
                    // Both ratebol and ratebur are modified by temperature (tfrt) and ratebur is also affected by water stress (wfdb).

                    // rate of burr growth, g per boll per day.
                    let ratebur = if AgeOfBoll[k][l][m] >= 22. {
                        0.
                    } else {
                        // Compute wfdb for the effect of water stress on burr growth rate.
                        // wfdb is the effect of water stress on rate of burr growth.
                        let wfdb = vpotfrt[0] + vpotfrt[1] * WaterStress;
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
                    PotGroBolls[k][l][m] = ratebol * FruitFraction[k][l][m];
                    PotGroBurrs[k][l][m] = ratebur * FruitFraction[k][l][m];
                    // Sum potential growth rates of bolls and burrs as PotGroAllBolls and PotGroAllBurrs, respectively.
                    PotGroAllBolls += PotGroBolls[k][l][m];
                    PotGroAllBurrs += PotGroBurrs[k][l][m];
                }
                // If these are not green bolls, their potential growth is 0. End loop.
                else {
                    PotGroBolls[k][l][m] = 0.;
                    PotGroBurrs[k][l][m] = 0.;
                }
            }
        }
    }
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
unsafe fn PotentialStemGrowth(stemnew: f64) -> f64 {
    // There are two periods for computation of potential stem growth:
    // (1) Before the appearance of a square on the third fruiting branch.
    // Potential stem growth is a functon of plant age (Kday, days from emergence).
    if FruitingCode[0][2][0] == 0 {
        VarPar[12] * (VarPar[13] + VarPar[14] * Kday as f64)
    }
    // (2) After the appearance of a square on the third fruiting branch.
    // It is assumed that all stem tissue that is more than 32 days old is not active.
    // Potential stem growth is a function of active stem tissue weight (stemnew), and plant density (denfac).
    else {
        // effect of plant density on stem growth rate.
        let denfac = 1. - VarPar[15] * (1. - DensityFactor);
        fmax(denfac, 0.2) * VarPar[16] * (VarPar[17] + VarPar[18] * stemnew)
    }
}
