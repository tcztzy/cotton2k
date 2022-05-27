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
    cgind, dl, maxk, maxl, nk, nl, pixcon, wk, ActualRootGrowth, CarbonAllocatedForRootGrowth,
    CultivationDate, CumPlantNLoss, DailyRootLoss, DayEmerge, Daynum, DepthLastRootLayer,
    ExtraCarbon, InitiateLateralRoots, LastTaprootLayer, LateralRootFlag, LateralRootGrowthLeft,
    LateralRootGrowthRight, NumLayersWithRoots, NumRootAgeGroups, PerPlantArea, PixInPlants,
    PlantRowColumn, PoreSpace, PotGroRoots, RedistRootNewGrowth, RootAge, RootAging,
    RootCultivation, RootDeath, RootGroFactor, RootImpedance, RootNConc, RootNitrogen,
    RootSummation, RootWeight, RootWeightLoss, RowSpace, SoilAirOnRootGrowth,
    SoilMechanicResistance, SoilNitrateOnRootGrowth, SoilPsi, SoilTemOnRootGrowth,
    SoilTempDailyAvrg, SoilWaterOnRootGrowth, State, TapRootGrowth, TapRootLength, VolNo3NContent,
    VolWaterContent,
};
use ndarray::prelude::*;
use ndarray::{Array, Array2, Ix2};

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

pub trait RootGrowth {
    unsafe fn compute_actual_root_growth(&mut self, sumpdr: f64);
}
impl RootGrowth for State {
    /// This function calculates the actual root growth rate.
    /// It is called from function PlantGrowth().
    /// It calls the following functions:  InitiateLateralRoots(), LateralRootGrowthLeft(), LateralRootGrowthRight(), RedistRootNewGrowth(), RootAging(), RootCultivation(), RootDeath(), RootSummation(), TapRootGrowth().
    ///
    /// The following global variables are referenced here:
    /// CarbonAllocatedForRootGrowth, cgind, CultivationDate, DayEmerge,
    /// Daynum, DepthLastRootLayer, dl, ExtraCarbon, LateralRootFlag,
    /// LastTaprootLayer, nk, nl, NumLayersWithRoots, NumRootAgeGroups,
    /// PerPlantArea, pixcon, PlantRowColumn, PotGroRoots, RootAge, RootNConc,
    /// RowSpace, TapRootLength, wk.
    ///
    /// The following global variables are set here:
    ///
    /// ActualRootGrowth, CumPlantNLoss, DailyRootLoss, PixInPlants, RootNitrogen, RootWeight, RootWeightLoss.
    ///
    /// The following argument is used:
    /// * `sumpdr` - Sum of potential root growth rate for the whole slab.
    unsafe fn compute_actual_root_growth(&mut self, sumpdr: f64) {
        //     The following constant parameters are used:
        //     The index for the relative partitioning of root mass produced by new
        //     growth to class i.
        const RootGrowthIndex: [f64; 3] = [1.0, 0.0, 0.0];
        const rtminc: f64 = 0.0000001; // the threshold ratio of root mass capable of growth
                                       // to soil cell volume (g/cm3); when this threshold is
                                       // reached, a
                                       // part of root growth in this cell may be extended to adjoining cells.
                                       //     Assign zero to pavail if this is the day of emergence.
        if Daynum <= DayEmerge {
            self.pavail = 0.;
        }
        let mut adwr1 = Array2::zeros([maxl as usize, maxk as usize]); // actual growth rate from roots existing in this
                                                                       // soil cell.
                                                                       //     Assign zero to the arrays of actual root growth rate.
        for l in 0..nl as usize {
            for k in 0..nk as usize {
                ActualRootGrowth[l][k] = 0.;
            }
        }
        //     The amount of carbon allocated for root growth is calculated from
        //  CarbonAllocatedForRootGrowth, converted to g dry matter per slab, and
        //  added to previously allocated carbon that has not been used for growth
        //  (pavail). if there is no potential root growth, this will be stored in
        //  pavail. Otherwise, zero is assigned to pavail.
        if sumpdr <= 0. {
            self.pavail += CarbonAllocatedForRootGrowth * 0.01 * RowSpace / PerPlantArea;
            return;
        }
        // actual growth factor (ratio of available C to potential growth).
        //     The ratio of available C to potential root growth (actgf) is
        //     calculated.
        //  pavail (if not zero) is used here, and zeroed after being used.
        let actgf =
            (self.pavail + CarbonAllocatedForRootGrowth * 0.01 * RowSpace / PerPlantArea) / sumpdr;
        self.pavail = 0.;
        //
        for l in 0..NumLayersWithRoots as usize {
            for k in 0..nk as usize {
                //     adwr1(l,k), is proportional to the potential growth rate of
                //     roots in this cell.
                if RootAge[l][k] > 0. {
                    adwr1.slice_mut(s![l, k]).fill(PotGroRoots[l][k] * actgf);
                }
            }
        }
        //     If extra carbon is available, it is assumed to be added to the
        //     taproot.
        if ExtraCarbon > 0. {
            // available carbon for taproot growth, in g dry matter per slab.
            //  ExtraCarbon is converted to availt (g dry matter per slab).
            let availt = ExtraCarbon * 0.01 * RowSpace / PerPlantArea;
            let mut sdl = TapRootLength - DepthLastRootLayer;
            // distance from the tip of the taproot, cm.
            let mut tpwt = Array::<f64, Ix2>::zeros([maxl as usize, 2]); // proportionality factors for allocating added dry matter among taproot soil cells.
            let mut sumwt = 0.; // sum of the tpwt array.
                                //     Extra Carbon (availt) is added to soil cells with roots in the
                                //     columns immediately to the
                                //  left and to the right of the location of the plant row.
            for l in (0..LastTaprootLayer as usize + 1).rev() {
                //     The weighting factors for allocating the carbon (tpwt) are
                //     proportional to the volume
                //  of each soil cell and its distance (sdl) from the tip of the
                //  taproot.
                sdl += dl[l];
                tpwt.slice_mut(s![l, 0])
                    .fill(sdl * dl[l] * wk[PlantRowColumn as usize]);
                tpwt.slice_mut(s![l, 1])
                    .fill(sdl * dl[l] * wk[PlantRowColumn as usize + 1]);
                sumwt += tpwt.slice(s![l, 0..2usize]).sum();
            }
            //     The proportional amount of mass is added to the mass of the last
            //     (inactive)
            //  root class in each soil cell.
            for l in 0..LastTaprootLayer as usize + 1 {
                RootWeight[l][PlantRowColumn as usize][NumRootAgeGroups as usize - 1] +=
                    availt * tpwt.slice(s![l, 0]).first().unwrap() / sumwt;
                RootWeight[l][PlantRowColumn as usize + 1][NumRootAgeGroups as usize - 1] +=
                    availt * tpwt.slice(s![l, 1]).first().unwrap() / sumwt;
            }
        }
        //     Check each cell if the ratio of root weight capable of growth to cell
        //     volume (rtconc)
        //  exceeds the threshold rtminc, and call RedistRootNewGrowth() for this
        //  cell. Otherwise, all new growth is contained in the same cell, and the
        //  actual growth in this cell, ActualRootGrowth(l,k) will be equal to
        //  adwr1(l,k).
        for l in 0..nl as usize {
            for k in 0..nk as usize {
                if RootAge[l][k] > 0. {
                    let mut rtconc = 0.; // ratio of root weight capable of growth to
                                         // cell volume.
                    for i in 0..NumRootAgeGroups as usize {
                        rtconc += RootWeight[l][k][i] * cgind[i];
                    }
                    rtconc /= dl[l] * wk[k];
                    if rtconc > rtminc {
                        RedistRootNewGrowth(l as i32, k as i32, adwr1[[l, k]]);
                    } else {
                        ActualRootGrowth[l][k] += adwr1[[l, k]];
                    }
                }
            }
        }
        //     The new actual growth ActualRootGrowth(l,k) in each cell is
        //     partitioned among the root
        //  classes in it in proportion to the parameters RootGrowthIndex(i), and
        //  the previous values of RootWeight(k,l,i), and added to
        //  RootWeight(k,l,i).
        let mut sumind = 0.; // sum of the growth index for all classes in a cell.
        for i in 0..NumRootAgeGroups as usize {
            sumind += RootGrowthIndex[i];
        }
        for l in 0..NumLayersWithRoots as usize {
            for k in 0..nk as usize {
                if RootAge[l][k] > 0. {
                    let mut sumgr = 0.; // sum of growth index multiplied by root
                                        // weight, for all classes in a cell.
                    for i in 0..NumRootAgeGroups as usize {
                        sumgr += RootGrowthIndex[i] * RootWeight[l][k][i];
                    }
                    for i in 0..NumRootAgeGroups as usize {
                        if sumgr > 0. {
                            RootWeight[l][k][i] +=
                                ActualRootGrowth[l][k] * RootGrowthIndex[i] * RootWeight[l][k][i]
                                    / sumgr;
                        } else {
                            RootWeight[l][k][i] +=
                                ActualRootGrowth[l][k] * RootGrowthIndex[i] / sumind;
                        }
                    }
                }
            }
        }
        //     Call function TapRootGrowth() for taproot elongation, if the taproot
        //  has not already reached the bottom of the slab.
        if LastTaprootLayer < nl - 1 || TapRootLength < DepthLastRootLayer {
            TapRootGrowth();
        }
        //     Call functions for growth of lateral roots
        InitiateLateralRoots();
        for l in 0..LastTaprootLayer as usize {
            if LateralRootFlag[l] == 2 {
                LateralRootGrowthLeft(l as i32);
                LateralRootGrowthRight(l as i32);
            }
        }
        //     Initialize DailyRootLoss (weight of sloughed roots) for this day.
        DailyRootLoss = 0.;
        for l in 0..NumLayersWithRoots as usize {
            for k in 0..nk as usize {
                //     Check RootAge to determine if this soil cell contains roots,
                //     and then compute root
                //  aging and root death by calling RootAging() and RootDeath() for
                //  each soil cell with roots.
                if RootAge[l][k] > 0. {
                    RootAging(l as i32, k as i32);
                    RootDeath(l as i32, k as i32);
                }
            }
        }
        //     Check if cultivation is executed in this day and call
        //     RootCultivation().
        for j in 0..5 {
            if CultivationDate[j] == Daynum {
                RootCultivation(j as i32);
            }
        }
        //     Convert DailyRootLoss to g per plant units and add it to
        //     RootWeightLoss.
        DailyRootLoss = DailyRootLoss * 100. * PerPlantArea / RowSpace;
        RootWeightLoss += DailyRootLoss;
        //     Adjust RootNitrogen (root N content) and PixInPlants (plant Pix
        //     content)
        //  for loss by death of roots.
        RootNitrogen -= DailyRootLoss * RootNConc;
        CumPlantNLoss += DailyRootLoss * RootNConc;
        PixInPlants -= DailyRootLoss * pixcon;
        //     Call function RootSummation().
        RootSummation();
    }
}
