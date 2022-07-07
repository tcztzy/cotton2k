use crate::plant::root::{PotentialRootGrowth, RootGrowth};
use crate::{
    nadj, pixdz, ActualFruitGrowth, ActualLeafGrowth, ActualStemGrowth, AdjAddHeightRate,
    AgeOfPreFruNode, AgeOfSite, AirTemp, CarbonStress, DayInc, DensityFactor, DryMatterBalance,
    FruitingCode, Kday, KdayAdjust, LeafAreaIndex, NStressVeg, NumAdjustDays, NumFruitBranches,
    NumPreFruNodes, PerPlantArea, PlantHeight, PotGroAllRoots, PotGroStem, PotentialFruitGrowth,
    PotentialLeafGrowth, PotentialStemGrowth, RowSpace, StemWeight, TotalLeafArea,
    TotalPetioleWeight, TotalStemWeight, VarPar, WaterStressStem,
};

use super::Plant;

pub trait PlantGrowth {
    unsafe fn growth(&mut self);
    unsafe fn plant_height_increment(&mut self, x: f64) -> f64;
}

impl PlantGrowth for Plant {
    /// This function simulates the potential and actual growth of cotton plants.
    /// It is called from [Profile::simulate_this_day()], and it calls the following functions:
    /// ActualFruitGrowth(), ActualLeafGrowth(), ActualRootGrowth(),
    /// AddPlantHeight(), DryMatterBalance(), PotentialFruitGrowth(),
    /// PotentialLeafGrowth(), PotentialRootGrowth(), PotentialStemGrowth(), TotalLeafWeight().
    ///
    /// The following global variables are referenced here:
    /// ActualStemGrowth, DayInc, FirstSquare, FruitingCode, Kday, pixdz,
    /// PerPlantArea, RowSpace, WaterStressStem.
    ///
    /// The following global variables are set here:
    /// LeafAreaIndex, PlantHeight, PotGroAllRoots, PotGroStem, StemWeight,
    /// TotalLeafArea, TotalLeafWeight, TotalPetioleWeight, TotalStemWeight.
    unsafe fn growth(&mut self) {
        //     Call PotentialLeafGrowth() to compute potential growth rate of
        //     leaves.
        PotentialLeafGrowth();
        //     If it is after first square, call PotentialFruitGrowth() to compute
        //     potential
        //  growth rate of squares and bolls.
        if FruitingCode[0][0][0] > 0 {
            PotentialFruitGrowth();
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
        let sumpdr = PotentialRootGrowth();
        // Total potential growth rate of roots is converted from g per slab (sumpdr) to g per plant (PotGroAllRoots).
        PotGroAllRoots = sumpdr * 100. * PerPlantArea / RowSpace;
        // Limit PotGroAllRoots to (maxrtgr*PerPlantArea) g per plant per day.
        let maxrtgr = 0.045; // maximum possible potential root growth, g dm-2 day-1.
        if PotGroAllRoots > maxrtgr * PerPlantArea {
            PotGroAllRoots = maxrtgr * PerPlantArea;
        }
        // Call DryMatterBalance() to compute carbon balance, allocation of carbon to plant parts, and carbon stress. DryMatterBalance() also computes and returns the values of the following arguments:
        // cdleaf is carbohydrate requirement for leaf growth, g per plant per day. cdpet is carbohydrate requirement for petiole growth, g per plant per day. cdroot is carbohydrate requirement for root growth, g per plant per day. cdstem is carbohydrate requirement for stem growth, g per plant per day.
        let mut cdstem = 0.;
        let mut cdleaf = 0.;
        let mut cdpet = 0.;
        let mut cdroot = 0.;
        DryMatterBalance(&mut cdstem, &mut cdleaf, &mut cdpet, &mut cdroot);
        // If it is after first square, call ActualFruitGrowth() to compute actual growth rate of squares and bolls.
        if FruitingCode[0][0][0] > 0 {
            ActualFruitGrowth();
        }
        // Initialize TotalLeafWeight. It is assumed that cotyledons fall off at time of first square. Also initialize TotalLeafArea and TotalPetioleWeight.
        TotalPetioleWeight = 0.;
        // Call ActualLeafGrowth to compute actual growth rate of leaves and compute leaf area index.
        ActualLeafGrowth();
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
        self.compute_actual_root_growth(sumpdr);
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
