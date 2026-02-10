use crate::general_functions::GetFromClim;
use crate::{
    ginp, nadj, pixcon, AbscisedFruitSites, AbscisedLeafWeight, AbscissionLag, AdjGreenBollAbsc,
    AdjSquareAbsc, AgeOfBoll, AgeOfPreFruNode, AgeOfSite, BloomWeightLoss, BollWeight, BurrNConc,
    BurrNitrogen, BurrWeight, BurrWeightGreenBolls, BurrWeightOpenBolls, CarbonStress,
    CottonWeightGreenBolls, CottonWeightOpenBolls, CumPlantNLoss, DayFirstDef, DayInc, Daynum,
    FruitFraction, FruitingCode, Gintot, GreenBollsLost, Kday, KdayAdjust, LeafAge, LeafArea,
    LeafAreaIndex, LeafAreaMainStem, LeafAreaNodes, LeafAreaPreFru, LeafNConc, LeafNitrogen,
    LeafWeightMainStem, LeafWeightNodes, LeafWeightPreFru, NitrogenStress, NumAbscisedLeaves,
    NumAdjustDays, NumFruitBranches, NumFruitSites, NumGreenBolls, NumNodes, NumOpenBolls,
    NumPreFruNodes, NumSheddingTags, NumSquares, NumVegBranches, PerPlantArea, PercentDefoliation,
    PetioleNConc, PetioleNitrogen, PetioleWeightMainStem, PetioleWeightNodes, PetioleWeightPreFru,
    PixInPlants, ReserveC, SeedNConc, SeedNitrogen, ShedByCarbonStress, ShedByNitrogenStress,
    ShedByWaterStress, SquareNConc, SquareNitrogen, SquareWeight, TotalLeafArea, TotalLeafWeight,
    TotalPetioleWeight, TotalSquareWeight, VarPar, WaterStress, CLIMATE_METRIC_TMAX,
};

fn drop_leaf_age(lai: f64) -> f64 {
    140.0 - lai
}

pub fn leaf_abscission() {
    unsafe {
        if LeafAreaIndex <= 0.0001 {
            return;
        }

        let droplf = drop_leaf_age(LeafAreaIndex);
        pre_fruit_leaf_abscission(droplf);

        for k in 0..NumVegBranches as usize {
            let nbrch = NumFruitBranches[k];
            for l in 0..nbrch as usize {
                main_stem_leaf_abscission(k, l, droplf);
            }
        }

        if DayFirstDef > 0 && Daynum >= DayFirstDef {
            defoliation_leaf_abscission();
        }

        if ReserveC > 0.0 {
            let resratio = 0.20;
            let resmax = resratio * TotalLeafWeight();
            if ReserveC > resmax {
                AbscisedLeafWeight += ReserveC - resmax;
                ReserveC = resmax;
            }
        }

        LeafAreaIndex = TotalLeafArea() / PerPlantArea;
        if LeafAreaIndex < 0.0001 {
            LeafAreaIndex = 0.0001;
        }
    }
}

fn pre_fruit_leaf_abscission(droplf: f64) {
    unsafe {
        for j in 0..NumPreFruNodes as usize {
            if crate::FirstSquare > 0 {
                AgeOfPreFruNode[j] += DayInc;
            }
            if AgeOfPreFruNode[j] >= droplf && LeafAreaPreFru[j] > 0.0 && LeafAreaIndex > 0.1 {
                LeafArea[NodeLayerPreFru(j)] -= LeafAreaPreFru[j];
                AbscisedLeafWeight += LeafWeightPreFru[j] + PetioleWeightPreFru[j];
                TotalPetioleWeight -= PetioleWeightPreFru[j];
                PixInPlants -= LeafWeightPreFru[j] * pixcon;
                LeafNitrogen -= LeafWeightPreFru[j] * LeafNConc;
                PetioleNitrogen -= PetioleWeightPreFru[j] * PetioleNConc;
                CumPlantNLoss +=
                    LeafWeightPreFru[j] * LeafNConc + PetioleWeightPreFru[j] * PetioleNConc;
                LeafAreaPreFru[j] = 0.0;
                LeafWeightPreFru[j] = 0.0;
                PetioleWeightPreFru[j] = 0.0;
                if DayFirstDef > 0 && Daynum > DayFirstDef {
                    NumAbscisedLeaves += 1;
                }
            }
        }
    }
}

fn main_stem_leaf_abscission(k: usize, l: usize, droplf: f64) {
    unsafe {
        if LeafAge[k][l][0] > droplf && LeafAreaMainStem[k][l] > 0.0 && LeafAreaIndex > 0.1 {
            AbscisedLeafWeight += LeafWeightMainStem[k][l] + PetioleWeightMainStem[k][l];
            TotalPetioleWeight -= PetioleWeightMainStem[k][l];
            PixInPlants -= LeafWeightMainStem[k][l] * pixcon;
            LeafNitrogen -= LeafWeightMainStem[k][l] * LeafNConc;
            PetioleNitrogen -= PetioleWeightMainStem[k][l] * PetioleNConc;
            CumPlantNLoss +=
                LeafWeightMainStem[k][l] * LeafNConc + PetioleWeightMainStem[k][l] * PetioleNConc;
            LeafArea[crate::NodeLayer[k][l] as usize] -= LeafAreaMainStem[k][l];
            LeafAreaMainStem[k][l] = 0.0;
            LeafWeightMainStem[k][l] = 0.0;
            PetioleWeightMainStem[k][l] = 0.0;
            if DayFirstDef > 0 && Daynum > DayFirstDef {
                NumAbscisedLeaves += 1;
            }
        }

        let nnid = NumNodes[k][l];
        for m in 0..nnid as usize {
            fruit_node_leaf_abscission(k, l, m, droplf);
        }
    }
}

fn fruit_node_leaf_abscission(k: usize, l: usize, m: usize, droplf: f64) {
    unsafe {
        if LeafAge[k][l][m] >= droplf && LeafAreaNodes[k][l][m] > 0.0 && LeafAreaIndex > 0.1 {
            AbscisedLeafWeight += LeafWeightNodes[k][l][m] + PetioleWeightNodes[k][l][m];
            TotalPetioleWeight -= PetioleWeightNodes[k][l][m];
            PixInPlants -= LeafWeightNodes[k][l][m] * pixcon;
            LeafNitrogen -= LeafWeightNodes[k][l][m] * LeafNConc;
            PetioleNitrogen -= PetioleWeightNodes[k][l][m] * PetioleNConc;
            CumPlantNLoss +=
                LeafWeightNodes[k][l][m] * LeafNConc + PetioleWeightNodes[k][l][m] * PetioleNConc;
            LeafArea[crate::NodeLayer[k][l] as usize] -= LeafAreaNodes[k][l][m];
            LeafAreaNodes[k][l][m] = 0.0;
            LeafWeightNodes[k][l][m] = 0.0;
            PetioleWeightNodes[k][l][m] = 0.0;
            if DayFirstDef > 0 && Daynum > DayFirstDef {
                NumAbscisedLeaves += 1;
            }
        }
    }
}

fn defoliation_leaf_abscission() {
    unsafe {
        if Daynum == DayFirstDef {
            for j in 0..NumPreFruNodes as usize {
                if LeafAreaPreFru[j] > 0.0 {
                    LeafArea[NodeLayerPreFru(j)] -= LeafAreaPreFru[j];
                    LeafAreaPreFru[j] = 0.0;
                    AbscisedLeafWeight += LeafWeightPreFru[j] + PetioleWeightPreFru[j];
                    TotalPetioleWeight -= PetioleWeightPreFru[j];
                    PixInPlants -= LeafWeightPreFru[j] * pixcon;
                    LeafNitrogen -= LeafWeightPreFru[j] * LeafNConc;
                    PetioleNitrogen -= PetioleWeightPreFru[j] * PetioleNConc;
                    CumPlantNLoss +=
                        LeafWeightPreFru[j] * LeafNConc + PetioleWeightPreFru[j] * PetioleNConc;
                    LeafWeightPreFru[j] = 0.0;
                    PetioleWeightPreFru[j] = 0.0;
                }
            }
        }

        if Daynum <= DayFirstDef {
            return;
        }

        let mut sort_by_age = [0.0; 450];
        let mut indexk = [0; 450];
        let mut indexl = [0; 450];
        let mut indexm = [0; 450];

        let mut lefcnt: usize = 0;
        for k in 0..NumVegBranches as usize {
            let nbrch = NumFruitBranches[k];
            for l in 0..nbrch as usize {
                if LeafWeightMainStem[k][l] > 0.0 {
                    sort_by_age[lefcnt] = AgeOfSite[k][l][0];
                    indexk[lefcnt] = k as i32;
                    indexl[lefcnt] = l as i32;
                    indexm[lefcnt] = 66;
                    lefcnt += 1;
                }
                let nnid = NumNodes[k][l];
                for m in 0..nnid as usize {
                    if LeafWeightNodes[k][l][m] > 0.0 {
                        sort_by_age[lefcnt] = AgeOfSite[k][l][m];
                        indexk[lefcnt] = k as i32;
                        indexl[lefcnt] = l as i32;
                        indexm[lefcnt] = m as i32;
                        lefcnt += 1;
                    }
                }
            }
        }

        sort_array(
            lefcnt,
            &mut sort_by_age,
            &mut indexk,
            &mut indexl,
            &mut indexm,
        );

        let mut num_leaves_to_shed = (lefcnt as f64 * PercentDefoliation / 100.0) as i32;
        for i in 0..lefcnt {
            if num_leaves_to_shed <= 0 {
                break;
            }
            if sort_by_age[i] > 0.0 {
                let k = indexk[i] as usize;
                let l = indexl[i] as usize;
                let m = indexm[i];
                let pixlos;
                if m == 66 {
                    AbscisedLeafWeight += LeafWeightMainStem[k][l] + PetioleWeightMainStem[k][l];
                    TotalPetioleWeight -= PetioleWeightMainStem[k][l];
                    pixlos = LeafWeightMainStem[k][l] * pixcon;
                    LeafNitrogen -= LeafWeightMainStem[k][l] * LeafNConc;
                    PetioleNitrogen -= PetioleWeightMainStem[k][l] * PetioleNConc;
                    CumPlantNLoss += LeafWeightMainStem[k][l] * LeafNConc
                        + PetioleWeightMainStem[k][l] * PetioleNConc;
                    LeafArea[crate::NodeLayer[k][l] as usize] -= LeafAreaMainStem[k][l];
                    LeafAreaMainStem[k][l] = 0.0;
                    LeafWeightMainStem[k][l] = 0.0;
                    PetioleWeightMainStem[k][l] = 0.0;
                } else {
                    let m = m as usize;
                    AbscisedLeafWeight += LeafWeightNodes[k][l][m] + PetioleWeightNodes[k][l][m];
                    TotalPetioleWeight -= PetioleWeightNodes[k][l][m];
                    pixlos = LeafWeightNodes[k][l][m] * pixcon;
                    LeafNitrogen -= LeafWeightNodes[k][l][m] * LeafNConc;
                    PetioleNitrogen -= PetioleWeightNodes[k][l][m] * PetioleNConc;
                    CumPlantNLoss += LeafWeightNodes[k][l][m] * LeafNConc
                        + PetioleWeightNodes[k][l][m] * PetioleNConc;
                    LeafArea[crate::NodeLayer[k][l] as usize] -= LeafAreaNodes[k][l][m];
                    LeafAreaNodes[k][l][m] = 0.0;
                    LeafWeightNodes[k][l][m] = 0.0;
                    PetioleWeightNodes[k][l][m] = 0.0;
                }
                PixInPlants -= pixlos;
                num_leaves_to_shed -= 1;
                NumAbscisedLeaves += 1;
            }
        }
    }
}

fn sort_array(
    size: usize,
    in_data: &mut [f64; 450],
    indexk: &mut [i32; 450],
    indexl: &mut [i32; 450],
    indexm: &mut [i32; 450],
) {
    let mut outk = [0; 450];
    let mut outl = [0; 450];
    let mut outm = [0; 450];
    let mut out_data = [0.0; 450];

    for i in 0..size {
        let value = in_data[i];
        let n0 = indexk[i];
        let n1 = indexl[i];
        let n2 = indexm[i];
        for j in 0..size {
            if value > out_data[j] {
                for k in (j + 1..size).rev() {
                    out_data[k] = out_data[k - 1];
                    outk[k] = outk[k - 1];
                    outl[k] = outl[k - 1];
                    outm[k] = outm[k - 1];
                }
                out_data[j] = value;
                outk[j] = n0;
                outl[j] = n1;
                outm[j] = n2;
                break;
            }
        }
    }

    for i in 0..size {
        indexk[i] = outk[i];
        indexl[i] = outl[i];
        indexm[i] = outm[i];
        in_data[i] = out_data[i];
    }
}

pub fn fruiting_sites_abscission() {
    unsafe {
        const VABSFR: [f64; 9] = [21.0, 0.42, 30.0, 0.05, 6.0, 2.25, 0.60, 5.0, 0.20];

        NumSheddingTags += 1;
        if NumSheddingTags > 1 {
            for lt in (1..NumSheddingTags as usize).rev() {
                let ltm1 = lt - 1;
                ShedByCarbonStress[lt] = ShedByCarbonStress[ltm1];
                ShedByNitrogenStress[lt] = ShedByNitrogenStress[ltm1];
                ShedByWaterStress[lt] = ShedByWaterStress[ltm1];
                AbscissionLag[lt] = AbscissionLag[ltm1];
            }
        }

        if CarbonStress < VarPar[43] {
            ShedByCarbonStress[0] = (VarPar[43] - CarbonStress) / VarPar[43];
        } else {
            ShedByCarbonStress[0] = 0.0;
        }
        if NitrogenStress < VABSFR[1] {
            ShedByNitrogenStress[0] = (VABSFR[1] - NitrogenStress) / VABSFR[1];
        } else {
            ShedByNitrogenStress[0] = 0.0;
        }
        if WaterStress < VarPar[44] {
            ShedByWaterStress[0] = (VarPar[44] - WaterStress) / VarPar[44];
        } else {
            ShedByWaterStress[0] = 0.0;
        }
        AbscissionLag[0] = 0.01;

        let tmax = GetFromClim(CLIMATE_METRIC_TMAX, Daynum);
        for lt in 0..NumSheddingTags as usize {
            AbscissionLag[lt] += DayInc.max(0.40);
            if tmax > VABSFR[2] {
                AbscissionLag[lt] += (tmax - VABSFR[2]) * VABSFR[3];
            }
        }

        let mut idecr = 0;
        for lt in 0..NumSheddingTags as usize {
            if AbscissionLag[lt] >= VABSFR[4] || lt >= 20 {
                let gin1 = if Gintot > 0.0 { Gintot } else { ginp };
                for k in 0..NumVegBranches as usize {
                    let nbrch = NumFruitBranches[k];
                    for l in 0..nbrch as usize {
                        let nnid = NumNodes[k][l];
                        for m in 0..nnid as usize {
                            if FruitingCode[k][l][m] == 1
                                || FruitingCode[k][l][m] == 2
                                || FruitingCode[k][l][m] == 7
                            {
                                let abscission_ratio = site_abscission_ratio(k, l, m, lt);
                                if abscission_ratio > 0.0 {
                                    if FruitingCode[k][l][m] == 1 {
                                        square_abscission(k, l, m, abscission_ratio);
                                    } else {
                                        boll_abscission(k, l, m, abscission_ratio, gin1);
                                    }
                                }
                            }
                        }
                    }
                }

                ShedByCarbonStress[lt] = 0.0;
                ShedByNitrogenStress[lt] = 0.0;
                ShedByWaterStress[lt] = 0.0;
                AbscissionLag[lt] = 0.0;
                idecr += 1;
            }
        }

        NumSheddingTags -= idecr;
        if Kday > KdayAdjust && Kday <= (KdayAdjust + NumAdjustDays) {
            adjust_abscission();
        }

        compute_site_numbers();
    }
}

fn site_abscission_ratio(k: usize, l: usize, m: usize, lt: usize) -> f64 {
    unsafe {
        const VABSC: [f64; 5] = [21.0, 2.25, 0.60, 5.0, 0.20];

        let mut pabs = 0.0;
        let mut shedt = 0.0;
        if FruitingCode[k][l][m] == 1 {
            if AgeOfSite[k][l][m] < VABSC[3] {
                pabs = 0.0;
            } else {
                let xsqage = AgeOfSite[k][l][m] - VABSC[3];
                if xsqage >= VABSC[0] {
                    pabs = VarPar[46];
                } else {
                    pabs = VarPar[46]
                        + (VarPar[45] - VarPar[46])
                            * ((VABSC[0] - xsqage) / VABSC[0]).powf(VABSC[1]);
                }
            }
            shedt = 1.0 - (1.0 - ShedByCarbonStress[lt]) * (1.0 - ShedByNitrogenStress[lt]);
        } else if FruitingCode[k][l][m] == 7 && AgeOfBoll[k][l][m] <= VarPar[47] {
            pabs = VarPar[48];
            shedt =
                1.0 - (1.0 - ShedByCarbonStress[lt]) * (1.0 - VABSC[2] * ShedByNitrogenStress[lt]);
        } else if AgeOfBoll[k][l][m] > VarPar[47] && AgeOfBoll[k][l][m] <= (VarPar[47] + VarPar[49])
        {
            pabs = VarPar[48]
                - (VarPar[48] - VarPar[50]) * (AgeOfBoll[k][l][m] - VarPar[47]) / VarPar[49];
            shedt = 1.0
                - (1.0 - ShedByCarbonStress[lt])
                    * (1.0 - VABSC[4] * ShedByNitrogenStress[lt])
                    * (1.0 - ShedByWaterStress[lt]);
        } else if AgeOfBoll[k][l][m] > (VarPar[47] + VarPar[49])
            && AgeOfBoll[k][l][m] <= (VarPar[47] + 2.0 * VarPar[49])
        {
            pabs = VarPar[50] / VarPar[49] * (VarPar[47] + 2.0 * VarPar[49] - AgeOfBoll[k][l][m]);
            shedt = ShedByWaterStress[lt];
        } else if AgeOfBoll[k][l][m] > (VarPar[47] + 2.0 * VarPar[49]) {
            pabs = 0.0;
        }

        let mut abscission_ratio = pabs * shedt * DayInc;
        if abscission_ratio > 1.0 {
            abscission_ratio = 1.0;
        }
        abscission_ratio
    }
}

fn square_abscission(k: usize, l: usize, m: usize, abscission_ratio: f64) {
    unsafe {
        let wtlos = SquareWeight[k][l][m] * abscission_ratio;
        SquareNitrogen -= wtlos * SquareNConc;
        CumPlantNLoss += wtlos * SquareNConc;
        SquareWeight[k][l][m] -= wtlos;
        BloomWeightLoss += wtlos;
        TotalSquareWeight -= wtlos;
        FruitFraction[k][l][m] *= 1.0 - abscission_ratio;

        if FruitFraction[k][l][m] <= 0.001 {
            FruitFraction[k][l][m] = 0.0;
            SquareNitrogen -= SquareWeight[k][l][m] * SquareNConc;
            CumPlantNLoss += SquareWeight[k][l][m] * SquareNConc;
            BloomWeightLoss += SquareWeight[k][l][m];
            TotalSquareWeight -= SquareWeight[k][l][m];
            SquareWeight[k][l][m] = 0.0;
            FruitingCode[k][l][m] = 5;
        }
    }
}

fn boll_abscission(k: usize, l: usize, m: usize, abscission_ratio: f64, gin1: f64) {
    unsafe {
        SeedNitrogen -= BollWeight[k][l][m] * abscission_ratio * (1.0 - gin1) * SeedNConc;
        BurrNitrogen -= BurrWeight[k][l][m] * abscission_ratio * BurrNConc;
        CumPlantNLoss += BollWeight[k][l][m] * abscission_ratio * (1.0 - gin1) * SeedNConc;
        CumPlantNLoss += BurrWeight[k][l][m] * abscission_ratio * BurrNConc;
        PixInPlants -= (BollWeight[k][l][m] + BurrWeight[k][l][m]) * abscission_ratio * pixcon;
        GreenBollsLost += (BollWeight[k][l][m] + BurrWeight[k][l][m]) * abscission_ratio;
        CottonWeightGreenBolls -= BollWeight[k][l][m] * abscission_ratio;
        BurrWeightGreenBolls -= BurrWeight[k][l][m] * abscission_ratio;
        BollWeight[k][l][m] -= BollWeight[k][l][m] * abscission_ratio;
        BurrWeight[k][l][m] -= BurrWeight[k][l][m] * abscission_ratio;
        FruitFraction[k][l][m] -= FruitFraction[k][l][m] * abscission_ratio;

        if FruitFraction[k][l][m] <= 0.001 {
            FruitingCode[k][l][m] = 4;
            SeedNitrogen -= BollWeight[k][l][m] * (1.0 - gin1) * SeedNConc;
            BurrNitrogen -= BurrWeight[k][l][m] * BurrNConc;
            CumPlantNLoss += BollWeight[k][l][m] * (1.0 - gin1) * SeedNConc;
            CumPlantNLoss += BurrWeight[k][l][m] * BurrNConc;
            PixInPlants -= (BollWeight[k][l][m] + BurrWeight[k][l][m]) * pixcon;
            FruitFraction[k][l][m] = 0.0;
            CottonWeightGreenBolls -= BollWeight[k][l][m];
            BurrWeightGreenBolls -= BurrWeight[k][l][m];
            GreenBollsLost += BollWeight[k][l][m] + BurrWeight[k][l][m];
            BollWeight[k][l][m] = 0.0;
            BurrWeight[k][l][m] = 0.0;
        }
    }
}

fn adjust_abscission() {
    unsafe {
        let mut jx = [0_i32; 2];
        let mut abscsq = 0.0;
        if nadj[3] && AdjSquareAbsc > 0.0 {
            jx[0] = 1;
            abscsq = AdjSquareAbsc;
        }

        let mut abscgb = 0.0;
        if nadj[4] && AdjGreenBollAbsc > 0.0 {
            jx[1] = 1;
            let mut gbolnum = 0.0;
            for k in 0..NumVegBranches as usize {
                let nbrch = NumFruitBranches[k];
                for l in 0..nbrch as usize {
                    let nnid = NumNodes[k][l];
                    for m in 0..nnid as usize {
                        if FruitingCode[k][l][m] == 7 {
                            gbolnum += FruitFraction[k][l][m];
                        }
                    }
                }
            }
            if gbolnum > 0.0 {
                abscgb = AdjGreenBollAbsc * NumGreenBolls / gbolnum;
            }
        }

        let mut gin1 = 0.0;
        if jx[1] == 1 {
            gin1 = if Gintot > 0.0 { Gintot } else { ginp };
        }

        for k in 0..NumVegBranches as usize {
            let nbrch = NumFruitBranches[k];
            for l in 0..nbrch as usize {
                let nnid = NumNodes[k][l];
                for m in 0..nnid as usize {
                    if jx[0] == 1 && FruitingCode[k][l][m] == 1 {
                        adjust_square_abscission(k, l, m, abscsq);
                    }
                    if jx[1] == 1 && FruitingCode[k][l][m] == 7 {
                        adjust_young_boll_abscission(k, l, m, abscgb, gin1);
                    }
                }
            }
        }
    }
}

fn adjust_square_abscission(k: usize, l: usize, m: usize, abscsq: f64) {
    unsafe {
        let mut wtlos = SquareWeight[k][l][m] * abscsq;
        SquareNitrogen -= wtlos * SquareNConc;
        CumPlantNLoss += wtlos * SquareNConc;
        SquareWeight[k][l][m] -= wtlos;
        BloomWeightLoss += wtlos;
        TotalSquareWeight -= wtlos;
        FruitFraction[k][l][m] *= 1.0 - abscsq;

        if FruitFraction[k][l][m] <= 0.001 {
            FruitFraction[k][l][m] = 0.0;
            SquareNitrogen -= SquareWeight[k][l][m] * SquareNConc;
            CumPlantNLoss += SquareWeight[k][l][m] * SquareNConc;
            BloomWeightLoss += SquareWeight[k][l][m];
            TotalSquareWeight -= SquareWeight[k][l][m];
            SquareWeight[k][l][m] = 0.0;
            FruitingCode[k][l][m] = 5;
        }

        if FruitFraction[k][l][m] > 1.0 {
            wtlos = SquareWeight[k][l][m] * (1.0 - 1.0 / FruitFraction[k][l][m]);
            FruitFraction[k][l][m] = 1.0;
            SquareNitrogen -= wtlos * SquareNConc;
            CumPlantNLoss += wtlos * SquareNConc;
            SquareWeight[k][l][m] -= wtlos;
            BloomWeightLoss += wtlos;
            TotalSquareWeight -= wtlos;
        }
    }
}

fn adjust_young_boll_abscission(k: usize, l: usize, m: usize, abscgb: f64, gin1: f64) {
    unsafe {
        SeedNitrogen -= BollWeight[k][l][m] * abscgb * (1.0 - gin1) * SeedNConc;
        BurrNitrogen -= BurrWeight[k][l][m] * abscgb * BurrNConc;
        CumPlantNLoss += BollWeight[k][l][m] * abscgb * (1.0 - gin1) * SeedNConc;
        CumPlantNLoss += BurrWeight[k][l][m] * abscgb * BurrNConc;
        let pixlos = (BollWeight[k][l][m] + BurrWeight[k][l][m]) * abscgb * pixcon;
        PixInPlants -= pixlos;
        GreenBollsLost += (BollWeight[k][l][m] + BurrWeight[k][l][m]) * abscgb;
        CottonWeightGreenBolls -= BollWeight[k][l][m] * abscgb;
        BurrWeightGreenBolls -= BurrWeight[k][l][m] * abscgb;
        BollWeight[k][l][m] *= 1.0 - abscgb;
        BurrWeight[k][l][m] *= 1.0 - abscgb;
        FruitFraction[k][l][m] *= 1.0 - abscgb;

        adjust_boll_abscission(k, l, m, 2, gin1);
    }
}

fn adjust_boll_abscission(k: usize, l: usize, m: usize, jx: i32, gin1: f64) {
    unsafe {
        if FruitFraction[k][l][m] <= 0.001 {
            SeedNitrogen -= BollWeight[k][l][m] * (1.0 - gin1) * SeedNConc;
            BurrNitrogen -= BurrWeight[k][l][m] * BurrNConc;
            CumPlantNLoss += BollWeight[k][l][m] * (1.0 - gin1) * SeedNConc;
            CumPlantNLoss += BurrWeight[k][l][m] * BurrNConc;
            PixInPlants -= (BollWeight[k][l][m] + BurrWeight[k][l][m]) * pixcon;
            GreenBollsLost += BollWeight[k][l][m] + BurrWeight[k][l][m];

            if jx == 2 {
                CottonWeightGreenBolls -= BollWeight[k][l][m];
                BurrWeightGreenBolls -= BurrWeight[k][l][m];
            } else if jx == 3 {
                CottonWeightOpenBolls -= BollWeight[k][l][m];
                BurrWeightOpenBolls -= BurrWeight[k][l][m];
            }

            FruitFraction[k][l][m] = 0.0;
            BollWeight[k][l][m] = 0.0;
            BurrWeight[k][l][m] = 0.0;
            FruitingCode[k][l][m] = 4;
        }

        if FruitFraction[k][l][m] > 1.0 {
            let bolwlos = BollWeight[k][l][m] * (1.0 - 1.0 / FruitFraction[k][l][m]);
            let burwlos = BurrWeight[k][l][m] * (1.0 - 1.0 / FruitFraction[k][l][m]);
            FruitFraction[k][l][m] = 1.0;
            SeedNitrogen -= bolwlos * (1.0 - gin1) * SeedNConc;
            BurrNitrogen -= burwlos * BurrNConc;
            CumPlantNLoss += bolwlos * (1.0 - gin1) * SeedNConc;
            CumPlantNLoss += burwlos * BurrNConc;
            PixInPlants -= (bolwlos + burwlos) * pixcon;
            GreenBollsLost += bolwlos + burwlos;
            BollWeight[k][l][m] -= bolwlos;
            BurrWeight[k][l][m] -= burwlos;

            if jx == 2 {
                CottonWeightGreenBolls -= bolwlos;
                BurrWeightGreenBolls -= burwlos;
            } else if jx == 3 {
                CottonWeightOpenBolls -= bolwlos;
                BurrWeightOpenBolls -= burwlos;
            }
        }
    }
}

fn compute_site_numbers() {
    unsafe {
        NumSquares = 0.0;
        NumGreenBolls = 0.0;
        NumOpenBolls = 0.0;
        for k in 0..NumVegBranches as usize {
            let nbrch = NumFruitBranches[k];
            for l in 0..nbrch as usize {
                let nnid = NumNodes[k][l];
                for m in 0..nnid as usize {
                    if FruitingCode[k][l][m] == 1 {
                        NumSquares += FruitFraction[k][l][m];
                    } else if FruitingCode[k][l][m] == 7 || FruitingCode[k][l][m] == 2 {
                        NumGreenBolls += FruitFraction[k][l][m];
                    } else if FruitingCode[k][l][m] == 3 {
                        NumOpenBolls += FruitFraction[k][l][m];
                    }
                }
            }
        }

        AbscisedFruitSites = NumFruitSites as f64 - NumSquares - NumGreenBolls - NumOpenBolls;
    }
}

#[inline]
fn NodeLayerPreFru(index: usize) -> usize {
    unsafe { crate::NodeLayerPreFru[index] as usize }
}
