use crate::{
    pixdz, ActualFruitGrowth, ActualLeafGrowth, ActualStemGrowth, AddPlantHeight, AirTemp,
    ComputeActualRootGrowth, DayInc, DensityFactor, DryMatterBalance, FruitingCode, Kday,
    LeafAreaIndex, PerPlantArea, PlantHeight, PotGroAllRoots, PotGroStem, PotentialFruitGrowth,
    PotentialLeafGrowth, PotentialRootGrowth, PotentialStemGrowth, RowSpace, StemWeight,
    TotalLeafArea, TotalPetioleWeight, TotalStemWeight, WaterStressStem,
};
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
pub unsafe fn PlantGrowth() {
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
    PlantHeight += AddPlantHeight(denf2);
    // Call ActualRootGrowth() to compute actual root growth.
    ComputeActualRootGrowth(sumpdr);
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
        let mut tfd = (AirTemp[ihr] - p1) / p2; // the hourly contribution to physiological age.
        if tfd < 0. {
            tfd = 0.;
        }
        if tfd > p3 {
            tfd = p3;
        }
        dayfd += tfd;
    }
    return dayfd / 24.;
}
