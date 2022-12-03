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
    cgind, dl, gh2oc, impede, inrim, maxk, maxl, ncurve, nk, nl, pixcon, rlat1, rlat2, tstbd, wk,
    ActualRootGrowth, BulkDensity, CarbonAllocatedForRootGrowth, CultivationDate, CumPlantNLoss,
    DailyRootLoss, DayEmerge, Daynum, DepthLastRootLayer, ExtraCarbon, InitiateLateralRoots,
    LastTaprootLayer, LateralRootFlag, NumLayersWithRoots, NumRootAgeGroups, PerPlantArea,
    PixInPlants, PlantRowColumn, PoreSpace, PotGroRoots, RootAge, RootColNumLeft, RootColNumRight,
    RootCultivation, RootGroFactor, RootImpede, RootNConc, RootNitrogen, RootSummation, RootWeight,
    RootWeightLoss, RowSpace, SoilAirOnRootGrowth, SoilHorizonNum, SoilMechanicResistance,
    SoilNitrateOnRootGrowth, SoilPsi, SoilTempDailyAvrg, SoilWaterOnRootGrowth, TapRootLength,
    VolNo3NContent, VolWaterContent,
};
use ndarray::prelude::*;
use ndarray::{Array, Array2, Ix2};

use super::Plant;

/// Calculates soil mechanical impedance to root growth, rtimpd(l,k), for all soil cells.
///
/// It is called from PotentialRootGrowth().
/// The impedance is a function of bulk density and water content in each soil soil cell.
/// No changes have been made in the original GOSSYM code.
///
/// The following global variables are referenced here:
/// BulkDensity, gh2oc, impede, inrim, ncurve, nk, nl,
/// SoilHorizonNum, tstbd, VolWaterContent.
///
/// The following global variables are set here:    RootImpede.
unsafe fn RootImpedance() {
    for l in 0..nl as usize {
        let j = SoilHorizonNum[l] as usize;
        let Bd = BulkDensity[j]; // bulk density for this layer
        let jj = tstbd.iter().position(|&x| Bd <= x[0]).unwrap();
        let j1 = if jj > inrim as usize - 1 {
            inrim as usize - 1
        } else {
            jj
        };
        let j0 = jj - 1;
        for k in 0..nk as usize {
            let Vh2o = VolWaterContent[l][k] / Bd;
            let ik = gh2oc.iter().position(|&x| Vh2o <= x).unwrap();
            let i1 = if ik > ncurve as usize - 1 {
                ncurve as usize - 1
            } else {
                ik
            };
            let i0 = ik - 1;
            if j1 == 0 {
                if i1 == 0 || Vh2o <= gh2oc[i1] {
                    RootImpede[l][k] = impede[j1][i1];
                } else {
                    RootImpede[l][k] = impede[j1][i0]
                        - (impede[j1][i0] - impede[j1][i1]) * (Vh2o - gh2oc[i0])
                            / (gh2oc[i1] - gh2oc[i0]);
                }
            } else {
                if i1 == 0 || Vh2o <= gh2oc[i1] {
                    RootImpede[l][k] = impede[j0][i1]
                        - (impede[j0][i1] - impede[j1][i1]) * (tstbd[j0][i1] - Bd)
                            / (tstbd[j0][i1] - tstbd[j1][i1]);
                } else {
                    let temp1 = impede[j0][i1]
                        - (impede[j0][i1] - impede[j1][i1]) * (tstbd[j0][i1] - Bd)
                            / (tstbd[j0][i1] - tstbd[j1][i1]);
                    let temp2 = impede[j0][i0]
                        - (impede[j0][i0] - impede[j1][i1]) * (tstbd[j0][i0] - Bd)
                            / (tstbd[j0][i0] - tstbd[j1][i0]);
                    RootImpede[l][k] =
                        temp2 + (temp1 - temp2) * (Vh2o - gh2oc[i0]) / (gh2oc[i1] - gh2oc[i0]);
                }
            }
        }
    }
    //
}

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

impl RootGrowth for Plant {
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
        // The following constant parameters are used:
        // The index for the relative partitioning of root mass produced by new growth to class i.
        const RootGrowthIndex: [f64; 3] = [1.0, 0.0, 0.0];
        // the threshold ratio of root mass capable of growth to soil cell volume (g/cm3); when this threshold is reached, a part of root growth in this cell may be extended to adjoining cells.
        const rtminc: f64 = 0.0000001;
        // Assign zero to pavail if this is the day of emergence.
        if Daynum <= DayEmerge {
            self.pavail = 0.;
        }
        let mut adwr1 = Array2::zeros([maxl as usize, maxk as usize]); // actual growth rate from roots existing in this soil cell.
                                                                       // Assign zero to the arrays of actual root growth rate.
        for l in 0..nl as usize {
            for k in 0..nk as usize {
                ActualRootGrowth[l][k] = 0.;
            }
        }
        // The amount of carbon allocated for root growth is calculated from
        // CarbonAllocatedForRootGrowth, converted to g dry matter per slab, and
        // added to previously allocated carbon that has not been used for growth
        // (pavail). If there is no potential root growth, this will be stored in
        // pavail. Otherwise, zero is assigned to pavail.
        if sumpdr <= 0. {
            self.pavail += CarbonAllocatedForRootGrowth * 0.01 * RowSpace / PerPlantArea;
            return;
        }
        // actual growth factor (ratio of available C to potential growth).
        // The ratio of available C to potential root growth (actgf) is calculated.
        // pavail (if not zero) is used here, and zeroed after being used.
        let actgf =
            (self.pavail + CarbonAllocatedForRootGrowth * 0.01 * RowSpace / PerPlantArea) / sumpdr;
        self.pavail = 0.;
        for l in 0..NumLayersWithRoots as usize {
            for k in 0..nk as usize {
                // adwr1(l,k), is proportional to the potential growth rate of roots in this cell.
                if RootAge[l][k] > 0. {
                    adwr1.slice_mut(s![l, k]).fill(PotGroRoots[l][k] * actgf);
                }
            }
        }
        // If extra carbon is available, it is assumed to be added to the taproot.
        if ExtraCarbon > 0. {
            // available carbon for taproot growth, in g dry matter per slab.
            //  ExtraCarbon is converted to availt (g dry matter per slab).
            let availt = ExtraCarbon * 0.01 * RowSpace / PerPlantArea;
            let mut sdl = TapRootLength - DepthLastRootLayer;
            // distance from the tip of the taproot, cm.
            let mut tpwt = Array::<f64, Ix2>::zeros([maxl as usize, 2]); // proportionality factors for allocating added dry matter among taproot soil cells.
            let mut sumwt = 0.; // sum of the tpwt array.
                                // Extra Carbon (availt) is added to soil cells with roots in the columns immediately to
                                // the left and to the right of the location of the plant row.
            for l in (0..LastTaprootLayer as usize + 1).rev() {
                // The weighting factors for allocating the carbon (tpwt) are proportional to the
                // volume of each soil cell and its distance (sdl) from the tip of the taproot.
                sdl += dl[l];
                tpwt.slice_mut(s![l, 0])
                    .fill(sdl * dl[l] * wk[PlantRowColumn as usize]);
                tpwt.slice_mut(s![l, 1])
                    .fill(sdl * dl[l] * wk[PlantRowColumn as usize + 1]);
                sumwt += tpwt.slice(s![l, 0..2usize]).sum();
            }
            // The proportional amount of mass is added to the mass of the last (inactive) root class in each soil cell.
            for l in 0..LastTaprootLayer as usize + 1 {
                RootWeight[l][PlantRowColumn as usize][NumRootAgeGroups as usize - 1] +=
                    availt * tpwt.slice(s![l, 0]).first().unwrap() / sumwt;
                RootWeight[l][PlantRowColumn as usize + 1][NumRootAgeGroups as usize - 1] +=
                    availt * tpwt.slice(s![l, 1]).first().unwrap() / sumwt;
            }
        }
        // Check each cell if the ratio of root weight capable of growth to cell volume (rtconc)
        // exceeds the threshold rtminc, and call RedistRootNewGrowth() for this cell.
        // Otherwise, all new growth is contained in the same cell, and the actual growth in this
        // cell, ActualRootGrowth(l,k) will be equal to adwr1(l,k).
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
                        RedistRootNewGrowth(l, k, adwr1[[l, k]]);
                    } else {
                        ActualRootGrowth[l][k] += adwr1[[l, k]];
                    }
                }
            }
        }
        // The new actual growth ActualRootGrowth(l,k) in each cell is partitioned among the root
        // classes in it in proportion to the parameters RootGrowthIndex(i), and the previous values
        // of RootWeight(k,l,i), and added to RootWeight(k,l,i).
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
        // Call function TapRootGrowth() for taproot elongation, if the taproot has not already
        // reached the bottom of the slab.
        if LastTaprootLayer < nl - 1 || TapRootLength < DepthLastRootLayer {
            TapRootGrowth();
        }
        // Call functions for growth of lateral roots
        InitiateLateralRoots();
        for l in 0..LastTaprootLayer as usize {
            if LateralRootFlag[l] == 2 {
                LateralRootGrowthLeft(l);
                LateralRootGrowthRight(l);
            }
        }
        // Initialize DailyRootLoss (weight of sloughed roots) for this day.
        DailyRootLoss = 0.;
        for l in 0..NumLayersWithRoots as usize {
            for k in 0..nk as usize {
                // Check RootAge to determine if this soil cell contains roots,
                // and then compute root aging and root death by calling RootAging() and RootDeath()
                // for each soil cell with roots.
                if RootAge[l][k] > 0. {
                    RootAging(l, k);
                    RootDeath(l, k);
                }
            }
        }
        // Check if cultivation is executed in this day and call RootCultivation().
        for j in 0..5 {
            if CultivationDate[j] == Daynum {
                RootCultivation(j as i32);
            }
        }
        // Convert DailyRootLoss to g per plant units and add it to RootWeightLoss.
        DailyRootLoss = DailyRootLoss * 100. * PerPlantArea / RowSpace;
        RootWeightLoss += DailyRootLoss;
        // Adjust RootNitrogen (root N content) and PixInPlants (plant Pix content)
        // for loss by death of roots.
        RootNitrogen -= DailyRootLoss * RootNConc;
        CumPlantNLoss += DailyRootLoss * RootNConc;
        PixInPlants -= DailyRootLoss * pixcon;
        // Call function RootSummation().
        RootSummation();
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
unsafe fn RedistRootNewGrowth(l: usize, k: usize, addwt: f64) {
    // The following constant parameters are used.
    // These are relative factors for root growth to adjoining cells, downwards, sideways, and upwards, respectively.
    // These factors are relative to the volume of the soil cell from which growth originates.
    const rgfdn: f64 = 900.;
    const rgfsd: f64 = 600.;
    const rgfup: f64 = 10.;
    // Set the number of layer above and below this layer, and the number of columns to the right and to the left of this column.
    // layer above and below layer l.
    let lp1 = if l == nl as usize - 1 { l } else { l + 1 };
    let lm1 = if l == 0 { l } else { l - 1 };
    // column to the left and to the right of column k.
    let kp1 = if k == nk as usize - 1 { k } else { k + 1 };
    let km1 = if k == 0 { k } else { k - 1 };
    // Compute proportionality factors (efac1, efacl, efacr, efacu, efacd) as the product of RootGroFactor and the geotropic factors in the respective soil cells.
    // Note that the geotropic factors are relative to the volume of the soil cell.
    // Compute the sum srwp of the proportionality factors.
    // product of RootGroFactor and geotropic factor for this cell.
    let efac1 = dl[l] * wk[k] * RootGroFactor[l][k];
    // as efac1 for the cell to the left of this cell.
    let efacl = rgfsd * RootGroFactor[l][km1];
    // as efac1 for the cell to the right of this cell.
    let efacr = rgfsd * RootGroFactor[l][kp1];
    // as efac1 for the cell above this cell.
    let efacu = rgfup * RootGroFactor[lm1][k];
    // as efac1 for the cell below this cell.
    let efacd = rgfdn * RootGroFactor[lp1][k];
    // sum of all efac values.
    let srwp = efac1 + efacl + efacr + efacu + efacd;
    // If srwp is very small, all the added weight will be in the
    // same soil soil cell, and execution of this function is ended.
    if srwp < 1e-10 {
        ActualRootGrowth[l][k] = addwt;
        return;
    }
    // Allocate the added dry matter to this and the adjoining soil cells in proportion to the EFAC factors.
    ActualRootGrowth[l][k] += addwt * efac1 / srwp;
    ActualRootGrowth[l][km1] += addwt * efacl / srwp;
    ActualRootGrowth[l][kp1] += addwt * efacr / srwp;
    ActualRootGrowth[lm1][k] += addwt * efacu / srwp;
    ActualRootGrowth[lp1][k] += addwt * efacd / srwp;
    // If roots are growing into new soil soil cells, initialize their RootAge to 0.01.
    if RootAge[l][km1] == 0. {
        RootAge[l][km1] = 0.01;
    }
    if RootAge[l][kp1] == 0. {
        RootAge[l][kp1] = 0.01;
    }
    if RootAge[lm1][k] == 0. {
        RootAge[lm1][k] = 0.01;
    }
    // If this new compartment is in a new layer with roots, also initialize its RootColNumLeft and RootColNumRight values.
    if RootAge[lp1][k] == 0. && efacd > 0. {
        RootAge[lp1][k] = 0.01;
        if RootColNumLeft[lp1] == 0 || k < RootColNumLeft[lp1] as usize {
            RootColNumLeft[lp1] = k as i32;
        }
        if RootColNumRight[lp1] == 0 || k > RootColNumRight[lp1] as usize {
            RootColNumRight[lp1] = k as i32;
        }
    }
    // If this is in the location of the taproot, and the roots reach a new soil layer,
    // update the taproot parameters TapRootLength, DepthLastRootLayer, and LastTaprootLayer.
    if k == PlantRowColumn as usize || k == PlantRowColumn as usize + 1 {
        if lp1 > LastTaprootLayer as usize && efacd > 0. {
            TapRootLength = DepthLastRootLayer + 0.01;
            DepthLastRootLayer += dl[lp1];
            LastTaprootLayer = lp1 as i32;
        }
    }
    // Update NumLayersWithRoots, if necessary, and the values of RootColNumLeft and RootColNumRight for this layer.
    if NumLayersWithRoots <= l as i32 && efacd > 0. {
        NumLayersWithRoots = l as i32 + 1;
    }
    if km1 < RootColNumLeft[l] as usize {
        RootColNumLeft[l] = km1 as i32;
    }
    if kp1 > RootColNumRight[l] as usize {
        RootColNumRight[l] = kp1 as i32;
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
unsafe fn RootAging(l: usize, k: usize) {
    //     The following constant parameters are used:
    const thtrn: [f64; 2] = [20.0, 40.0]; // the time threshold, from the initial
                                          // penetration of roots to a soil cell, after which some of the root
                                          // mass of class i may be transferred into the next class (i+1).
    const trn: [f64; 2] = [0.0060, 0.0050]; //  the daily proportion of this transfer.
                                            //
                                            // daily average soil temperature (c) of soil cell.
    let stday = SoilTempDailyAvrg[l][k] - 273.161;
    RootAge[l][k] += SoilTemOnRootGrowth(stday);
    //
    for i in 0..2 {
        if RootAge[l][k] > thtrn[i] {
            // root mass transferred from one class to the next.
            let xtr = trn[i] * RootWeight[l][k][i];
            RootWeight[l][k][i + 1] += xtr;
            RootWeight[l][k][i] -= xtr;
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
unsafe fn RootDeath(l: usize, k: usize) {
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
        if RootAge[l][k] > thdth[i] {
            // the computed proportion of roots dying in each
            // class.
            let mut dthfac = dth[i];
            if VolWaterContent[l][k] >= PoreSpace[l] {
                dthfac = dthmax;
            } else {
                if i <= 1 && SoilPsi[l][k] <= psi0 {
                    dthfac += aa * (psi0 - SoilPsi[l][k]);
                }
                if dthfac > dthmax {
                    dthfac = dthmax;
                }
            }
            DailyRootLoss += RootWeight[l][k][i] * dthfac;
            RootWeight[l][k][i] -= RootWeight[l][k][i] * dthfac;
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
unsafe fn TapRootGrowth() {
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
    let klocp1 = PlantRowColumn as usize + 1;
    if VolWaterContent[LastTaprootLayer as usize][PlantRowColumn as usize]
        >= PoreSpace[LastTaprootLayer as usize]
        || VolWaterContent[LastTaprootLayer as usize][klocp1]
            >= PoreSpace[LastTaprootLayer as usize]
    {
        return;
    }
    //     Average soil resistance (avres) is computed at the root tip.
    // avres = average value of RootGroFactor for the two soil cells at the tip
    // of the taproot.
    let avres = (RootGroFactor[LastTaprootLayer as usize][PlantRowColumn as usize]
        + RootGroFactor[LastTaprootLayer as usize][klocp1])
        / 2.;
    //     It is assumed that a linear empirical function of avres controls the
    //     rate of
    //  taproot elongation. The potential elongation rate of the taproot is also
    //  modified by soil temperature (SoilTemOnRootGrowth function), soil
    //  resistance, and soil moisture near the root tip.
    //     Actual growth is added to the taproot length TapRootLength.
    // daily average soil temperature (C) at root tip.
    let stday = 0.5
        * (SoilTempDailyAvrg[LastTaprootLayer as usize][PlantRowColumn as usize]
            + SoilTempDailyAvrg[LastTaprootLayer as usize][klocp1])
        - 273.161;
    // added taproot length, cm
    let addtaprt = rtapr * (1. - p1 + avres * p1) * SoilTemOnRootGrowth(stday);
    TapRootLength += addtaprt;
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
    if LastTaprootLayer > nl - 2 || TapRootLength <= DepthLastRootLayer {
        return;
    }
    //     The following is executed when the taproot reaches a new soil layer.
    LastTaprootLayer += 1;
    DepthLastRootLayer += dl[LastTaprootLayer as usize];
    if LastTaprootLayer > NumLayersWithRoots - 1 {
        NumLayersWithRoots = LastTaprootLayer + 1;
        if NumLayersWithRoots > nl {
            NumLayersWithRoots = nl;
        }
    }
    if RootColNumLeft[LastTaprootLayer as usize] == 0
        || RootColNumLeft[LastTaprootLayer as usize] > PlantRowColumn
    {
        RootColNumLeft[LastTaprootLayer as usize] = PlantRowColumn;
    }
    if RootColNumRight[LastTaprootLayer as usize] == 0
        || RootColNumRight[LastTaprootLayer as usize] < klocp1 as i32
    {
        RootColNumRight[LastTaprootLayer as usize] = klocp1 as i32;
    }
    //     RootAge is initialized for these soil cells.
    RootAge[LastTaprootLayer as usize][PlantRowColumn as usize] = 0.01;
    RootAge[LastTaprootLayer as usize][klocp1] = 0.01;
    //     Some of the mass of class 1 roots is transferred downwards to
    //  the new cells. The transferred mass is proportional to 2 cm of
    //  layer width, but it is not more than half the existing mass in the
    //  last layer.
    for i in 0..NumRootAgeGroups as usize {
        // root mass transferred to the cell below when the
        // elongating taproot reaches a new soil layer.
        // first column
        let mut tran = RootWeight[LastTaprootLayer as usize - 1][PlantRowColumn as usize][i] * 2.
            / dl[LastTaprootLayer as usize - 1];
        if tran > 0.5 * RootWeight[LastTaprootLayer as usize - 1][PlantRowColumn as usize][i] {
            tran = 0.5 * RootWeight[LastTaprootLayer as usize - 1][PlantRowColumn as usize][i];
        }
        RootWeight[LastTaprootLayer as usize][PlantRowColumn as usize][i] += tran;
        RootWeight[LastTaprootLayer as usize - 1][PlantRowColumn as usize][i] -= tran;
        // second column
        tran = RootWeight[LastTaprootLayer as usize - 1][klocp1][i] * 2.
            / dl[LastTaprootLayer as usize - 1];
        if tran > 0.5 * RootWeight[LastTaprootLayer as usize - 1][klocp1][i] {
            tran = 0.5 * RootWeight[LastTaprootLayer as usize - 1][klocp1][i];
        }
        RootWeight[LastTaprootLayer as usize][klocp1][i] += tran;
        RootWeight[LastTaprootLayer as usize - 1][klocp1][i] -= tran;
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
unsafe fn LateralRootGrowthLeft(l: usize) {
    //     The following constant parameters are used:
    const p1: f64 = 0.10; // constant parameter.
    const rlatr: f64 = 3.6; // potential growth rate of lateral roots, cm/day.
    const rtran: f64 = 0.2; // the ratio of root mass transferred to a new soil
                            // soil cell, when a lateral root grows into it.
                            //     On its initiation, lateral root length is assumed to be equal to the
                            //  width of a soil column soil cell at the location of the taproot.
    if rlat1[l] <= 0. {
        rlat1[l] = wk[PlantRowColumn as usize];
    }
    // daily average soil temperature (C) at root tip.
    let stday = SoilTempDailyAvrg[l][PlantRowColumn as usize] - 273.161;
    // the effect of soil temperature on root growth.
    let temprg = SoilTemOnRootGrowth(stday);
    //     Define the column with the tip of this lateral root (ktip)
    let mut ktip = 0usize; // column with the tips of the laterals to the left
    let mut sumwk = 0.; // summation of columns width
    for k in (0..PlantRowColumn as usize + 1).rev() {
        sumwk += wk[k];
        if sumwk >= rlat1[l] {
            ktip = k;
            break;
        }
    }
    //     Compute growth of the lateral root to the left. Potential
    //  growth rate (u) is modified by the soil temperature function,
    //  and the linearly modified effect of soil resistance (RootGroFactor).
    //     Lateral root elongation does not occur in water logged soil.
    if VolWaterContent[l][ktip] < PoreSpace[l] {
        rlat1[l] += rlatr * temprg * (1. - p1 + RootGroFactor[l][ktip] * p1);
        //     If the lateral reaches a new soil soil cell: a proportion (tran)
        //     of
        //	mass of roots is transferred to the new soil cell.
        if rlat1[l] > sumwk && ktip > 0 {
            // column into which the tip of the lateral grows to
            // left.
            let newktip = ktip - 1;
            for i in 0..NumRootAgeGroups as usize {
                let tran = RootWeight[l][ktip][i] * rtran;
                RootWeight[l][ktip][i] -= tran;
                RootWeight[l][newktip][i] += tran;
            }
            //     RootAge is initialized for this soil cell.
            //     RootColNumLeft of this layer is redefined.
            if RootAge[l][newktip] == 0. {
                RootAge[l][newktip] = 0.01;
            }
            if newktip < RootColNumLeft[l] as usize {
                RootColNumLeft[l] = newktip as i32;
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
unsafe fn LateralRootGrowthRight(l: usize) {
    //     The following constant parameters are used:
    const p1: f64 = 0.10; // constant parameter.
    const rlatr: f64 = 3.6; // potential growth rate of lateral roots, cm/day.
    const rtran: f64 = 0.2; // the ratio of root mass transferred to a new soil
                            // soil cell, when a lateral root grows into it.
                            //     On its initiation, lateral root length is assumed to be equal to the
                            //     width
                            //  of a soil column soil cell at the location of the taproot.
    let klocp1 = PlantRowColumn as usize + 1;
    if rlat2[l] <= 0. {
        rlat2[l] = wk[klocp1];
    }
    // daily average soil temperature (C) at root tip.
    let stday = SoilTempDailyAvrg[l][klocp1] - 273.161;
    // the effect of soil temperature on root growth.
    let temprg = SoilTemOnRootGrowth(stday);
    // define the column with the tip of this lateral root (ktip)
    let mut ktip = 0usize; // column with the tips of the laterals to the right
    let mut sumwk = 0.;
    for k in klocp1..nk as usize {
        sumwk += wk[k];
        if sumwk >= rlat2[l] {
            ktip = k;
            break;
        }
    }
    //     Compute growth of the lateral root to the right. Potential
    //  growth rate is modified by the soil temperature function,
    //  and the linearly modified effect of soil resistance (RootGroFactor).
    //     Lateral root elongation does not occur in water logged soil.
    if VolWaterContent[l][ktip] < PoreSpace[l] {
        rlat2[l] += rlatr * temprg * (1. - p1 + RootGroFactor[l][ktip] * p1);
        //     If the lateral reaches a new soil soil cell: a proportion (tran)
        //     of
        //	mass of roots is transferred to the new soil cell.
        if rlat2[l] > sumwk && ktip < nk as usize - 1 {
            // column into which the tip of the lateral grows to left.
            let newktip = ktip + 1; // column into which the tip of the lateral grows to left.
            for i in 0..NumRootAgeGroups as usize {
                let tran = RootWeight[l][ktip][i] * rtran;
                RootWeight[l][ktip][i] -= tran;
                RootWeight[l][newktip][i] += tran;
            }
            //     RootAge is initialized for this soil cell.
            //     RootColNumLeft of this layer is redefined.
            if RootAge[l][newktip] == 0. {
                RootAge[l][newktip] = 0.01;
            }
            if newktip > RootColNumRight[l] as usize {
                RootColNumRight[l] = newktip as i32;
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
unsafe fn SoilTemOnRootGrowth(t: f64) -> f64 {
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
