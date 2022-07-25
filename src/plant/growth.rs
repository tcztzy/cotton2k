use crate::atmosphere::Atmosphere;
use crate::plant::root::{PotentialRootGrowth, RootGrowth};
use crate::utils::fmax;
use crate::{
    nadj, pixdz, ActualFruitGrowth, ActualLeafGrowth, ActualStemGrowth, AdjAddHeightRate,
    AgeOfBoll, AgeOfPreFruNode, AgeOfSite, AirTemp, CarbonStress, DayInc, DayTimeTemp,
    DensityFactor, DryMatterBalance, FruitFraction, FruitingCode, Kday, KdayAdjust, LeafAreaIndex,
    NStressVeg, NightTimeTemp, NumAdjustDays, NumFruitBranches, NumNodes, NumPreFruNodes,
    NumVegBranches, PerPlantArea, PlantHeight, PotGroAllBolls, PotGroAllBurrs, PotGroAllRoots,
    PotGroAllSquares, PotGroBolls, PotGroBurrs, PotGroSquares, PotGroStem, PotentialLeafGrowth,
    RowSpace, StemWeight, TotalLeafArea, TotalPetioleWeight, TotalStemWeight, VarPar, WaterStress,
    WaterStressStem,
};

use super::Plant;
use crate::atmosphere::num_hours;
use chrono::Duration;

pub trait PlantGrowth {
    unsafe fn grow(&mut self, atmosphere: Atmosphere);
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
    unsafe fn grow(&mut self, atmosphere: Atmosphere) {
        //     Call PotentialLeafGrowth() to compute potential growth rate of
        //     leaves.
        PotentialLeafGrowth();
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
