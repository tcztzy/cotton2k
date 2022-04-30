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
//  at the start of the simulation (see their code in file gettinginput_2.cpp) .
//     PotentialRootGrowth() and ActualRootGrowth() are called each day from
//     PlantGrowth(). PotentialRootGrowth() calls RootImpedance(),
//     SoilMechanicResistance(), SoilAirOnRootGrowth(),
//  SoilNitrateOnRootGrowth(), SoilTemOnRootGrowth(), SoilWaterOnRootGrowth().
//     ActualRootGrowth() calls RedistRootNewGrowth(), TapRootGrowth(),
//     LateralRootGrowth(),
//  RootAging(), RootDeath(), RootCultivation(), RootSummation().
use crate::{
    cgind, nk, NumLayersWithRoots, NumRootAgeGroups, PerPlantArea, PoreSpace, PotGroRoots, RootAge,
    RootGroFactor, RootImpedance, RootWeight, SoilAirOnRootGrowth, SoilMechanicResistance,
    SoilNitrateOnRootGrowth, SoilPsi, SoilTemOnRootGrowth, SoilTempDailyAvrg,
    SoilWaterOnRootGrowth, VolNo3NContent, VolWaterContent,
};

/// Calculates the potential root growth rate.
/// The return value is the sum of potential root growth rates for the whole slab (sumpdr).
/// It is called from PlantGrowth(). It calls: RootImpedance(),
/// SoilNitrateOnRootGrowth(), SoilAirOnRootGrowth(), SoilMechanicResistance(),
/// SoilTemOnRootGrowth() and SoilWaterOnRootGrowth().
///
/// The following global variables are referenced here:
///
/// cgind, Date, Daynum, NumLayersWithRoots, NumRootAgeGroups, nk,
/// PerPlantArea, PoreSpace, RootAge, RootWeight. SoilPsi,
/// SoilTempDailyAvrg, VolNo3NContent, VolWaterContent.
///
/// The following global variables are set here:    PotGroRoots, RootGroFactor
pub unsafe fn PotentialRootGrowth() -> f64 {
    // The following constant parameter is used:
    // potential relative growth rate of the roots (g/g/day).
    const rgfac: f64 = 0.36;
    // Initialize to zero the PotGroRoots array.
    for l in 0..NumLayersWithRoots as usize {
        for k in 0..nk as usize {
            PotGroRoots[l][k] = 0.;
        }
    }
    RootImpedance();
    let mut sumpdr = 0.; // sum of potential root growth rate for the whole slab
    for l in 0..NumLayersWithRoots as usize {
        for k in 0..nk as usize {
            // Check if this soil cell contains roots (if RootAge is greater than 0), and execute the following if this is true.
            // In each soil cell with roots, the root weight capable of growth rtwtcg is computed as the sum of RootWeight[l][k][i] * cgind[i] for all root classes.
            if RootAge[l][k] > 0. {
                let mut rtwtcg = 0.; // root weight capable of growth in a soil soil cell.
                for i in 0..NumRootAgeGroups as usize {
                    rtwtcg += RootWeight[l][k][i] * cgind[i];
                }
                // Compute the temperature factor for root growth by calling function SoilTemOnRootGrowth() for this layer.
                // soil temperature, C, this day's average for this cell.
                let stday = SoilTempDailyAvrg[l][k] - 273.161;
                // effect of soil temperature on root growth.
                let temprg = SoilTemOnRootGrowth(stday);
                // Compute soil mechanical resistance for each soil cell by calling SoilMechanicResistance{},
                // the effect of soil aeration on root growth by calling SoilAirOnRootGrowth(),
                // and the effect of soil nitrate on root growth by calling SoilNitrateOnRootGrowth().
                // effect of soil mechanical resistance on root growth (returned from SoilMechanicResistance).
                let rtpct = SoilMechanicResistance(l as i32, k as i32);
                // effect of oxygen deficiency on root growth (returned from SoilAirOnRootGrowth).
                let rtrdo = SoilAirOnRootGrowth(SoilPsi[l][k], PoreSpace[l], VolWaterContent[l][k]);
                // effect of nitrate deficiency on root growth (returned from SoilNitrateOnRootGrowth).
                let rtrdn = SoilNitrateOnRootGrowth(VolNo3NContent[l][k]);
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
                let rtpsi = SoilWaterOnRootGrowth(SoilPsi[l][k]);
                RootGroFactor[l][k] = rtpsi * minres;
                PotGroRoots[l][k] =
                    rtwtcg * rgfac * temprg * RootGroFactor[l][k] * PerPlantArea / 19.6;
                sumpdr += PotGroRoots[l][k];
            }
        }
    }
    return sumpdr;
}
