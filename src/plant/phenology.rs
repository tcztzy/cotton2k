use crate::general_functions::GetFromClim;
use crate::plant::abscission::{fruiting_sites_abscission, leaf_abscission};
use crate::{
    bPollinSwitch, ginp, nadj, pixcon, pixdn, AbscisedLeafWeight, AdjAddMSNodesRate,
    AdjAddSitesRate, AgeOfBoll, AgeOfPreFruNode, AgeOfSite, AvrgDailyTemp, AvrgNodeTemper,
    BloomWeightLoss, BollWeight, BurrNitrogen, BurrWeight, BurrWeightGreenBolls,
    BurrWeightOpenBolls, CarbonStress, CottonWeightGreenBolls, CottonWeightOpenBolls,
    CumPlantNLoss, DayEmerge, DayFirstDef, DayInc, Daynum, DelayNewFruBranch, DelayNewNode,
    DensityFactor, FirstBloom, FirstSquare, FruitFraction, FruitGrowthRatio, FruitingCode, Gintot,
    Kday, KdayAdjust, LeafAge, LeafAreaIndex, LeafAreaMainStem, LeafAreaNodes, LeafAreaPreFru,
    LeafNitrogen, LeafWeightAreaRatio, LeafWeightMainStem, LeafWeightNodes, LeafWeightPreFru,
    LintYield, NStressFruiting, NStressVeg, NodeLayer, NodeLayerPreFru, NumAdjustDays,
    NumFruitBranches, NumFruitSites, NumNodes, NumOpenBolls, NumPreFruNodes, NumVegBranches,
    PerPlantArea, PixInPlants, PlantHeight, PlantPopulation, SeedNitrogen, SquareNConc,
    SquareNitrogen, SquareWeight, StemNitrogen, TotalLeafWeight, TotalSquareWeight,
    TotalStemWeight, VarPar, WaterStress, CLIMATE_METRIC_TMAX, CLIMATE_METRIC_TMIN,
};

static mut FIB_LENGTH: f64 = 0.0;
static mut FIB_STRENGTH: f64 = 0.0;
static mut PHEN_DELAY_BY_N_STRESS: f64 = 0.0;
static mut NWFL: i32 = 0;
static mut DAYS_TO_1ST_SQUARE: f64 = 0.0;
static mut BOLTMP: [[[f64; 5]; 30]; 3] = [[[0.0; 5]; 30]; 3];
static mut DAYS_TO_1ST_AVTEMP: f64 = 0.0;
static mut DAYS_TO_1ST_SUMSTRS: f64 = 0.0;

pub unsafe fn cotton_phenology() {
    const VPHENO: [f64; 8] = [0.65, -0.83, -1.67, -0.25, -0.75, 10.0, 15.0, 7.10];

    NumFruitSites = 0;
    let stem_n_ratio = StemNitrogen / TotalStemWeight;

    PHEN_DELAY_BY_N_STRESS = VPHENO[0] * (1.0 - NStressVeg);
    PHEN_DELAY_BY_N_STRESS = PHEN_DELAY_BY_N_STRESS.clamp(0.0, 1.0);

    let mut delay_veg_by_c_stress =
        VarPar[27] + CarbonStress * (VPHENO[3] + VPHENO[4] * CarbonStress);
    delay_veg_by_c_stress = delay_veg_by_c_stress.clamp(0.0, 1.0);

    let mut delay_frt_by_c_stress =
        VarPar[28] + CarbonStress * (VPHENO[1] + VPHENO[2] * CarbonStress);
    if delay_frt_by_c_stress > VarPar[29] {
        delay_frt_by_c_stress = VarPar[29];
    }
    if delay_frt_by_c_stress < 0.0 {
        delay_frt_by_c_stress = 0.0;
    }

    if FirstSquare <= 0 {
        DAYS_TO_1ST_SQUARE = days_to_first_square();
        pre_fruiting_node(stem_n_ratio);
        if Kday >= DAYS_TO_1ST_SQUARE as i32 {
            FirstSquare = Daynum;
            create_first_square(stem_n_ratio);
        } else {
            leaf_abscission();
            return;
        }
    }

    if NumVegBranches == 1 && PerPlantArea >= VPHENO[5] {
        add_vegetative_branch(delay_veg_by_c_stress, stem_n_ratio, DAYS_TO_1ST_SQUARE);
    }
    if NumVegBranches == 2 && PerPlantArea >= VPHENO[6] {
        add_vegetative_branch(delay_veg_by_c_stress, stem_n_ratio, DAYS_TO_1ST_SQUARE);
    }

    let mut nidmax = (VPHENO[7] * DensityFactor + 0.5) as i32;
    if nidmax > 5 {
        nidmax = 5;
    }

    let mut nwfl = NWFL;
    for k in 0..NumVegBranches as usize {
        if NumFruitBranches[k] < 30 {
            add_fruiting_branch(k, delay_veg_by_c_stress, stem_n_ratio);
        }

        for l in 0..NumFruitBranches[k] as usize {
            if NumNodes[k][l] < nidmax {
                add_fruiting_node(k, l, delay_frt_by_c_stress, stem_n_ratio);
            }

            for m in 0..NumNodes[k][l] as usize {
                fruiting_site(k, l, m, &mut nwfl);
            }
        }
    }
    NWFL = nwfl;

    fruiting_sites_abscission();
    leaf_abscission();
}

unsafe fn pre_fruiting_node(stem_n_ratio: f64) {
    const MAX_AGE_PRE_FR_NODE: f64 = 66.0;

    let last_index = (NumPreFruNodes - 1) as usize;
    if AgeOfPreFruNode[last_index] > MAX_AGE_PRE_FR_NODE {
        return;
    }

    for j in 0..NumPreFruNodes as usize {
        AgeOfPreFruNode[j] += DayInc;
    }

    if NumPreFruNodes >= 9 {
        return;
    }

    let mut time_to_next_pre_fru_node = VarPar[31];
    if NumPreFruNodes <= 2 {
        time_to_next_pre_fru_node *= VarPar[32];
    } else if NumPreFruNodes == 3 {
        time_to_next_pre_fru_node *= VarPar[33];
    }

    let last_index = (NumPreFruNodes - 1) as usize;
    if AgeOfPreFruNode[last_index] >= time_to_next_pre_fru_node {
        NumPreFruNodes += 1;
        let new_index = (NumPreFruNodes - 1) as usize;
        NodeLayerPreFru[new_index] = ((PlantHeight / 5.0) as i32).clamp(0, 19);
        LeafAreaPreFru[new_index] = VarPar[34];
        LeafWeightPreFru[new_index] = LeafAreaPreFru[new_index] * LeafWeightAreaRatio;
        TotalStemWeight -= LeafWeightPreFru[new_index];
        LeafNitrogen += LeafWeightPreFru[new_index] * stem_n_ratio;
        StemNitrogen -= LeafWeightPreFru[new_index] * stem_n_ratio;
    }
}

unsafe fn days_to_first_square() -> f64 {
    let p1 = 34.0;
    let p2 = 132.2;
    let p3 = -7.0;
    let p4 = 0.125;
    let p5 = 0.08;
    let p6 = 0.30;

    if Daynum <= DayEmerge {
        DAYS_TO_1ST_AVTEMP = AvrgDailyTemp;
        DAYS_TO_1ST_SUMSTRS = 0.0;
    }

    DAYS_TO_1ST_AVTEMP = ((Kday - 1) as f64 * DAYS_TO_1ST_AVTEMP + AvrgDailyTemp) / Kday as f64;
    if DAYS_TO_1ST_AVTEMP > p1 {
        DAYS_TO_1ST_AVTEMP = p1;
    }

    let mut tsq1 = p2 + DAYS_TO_1ST_AVTEMP * (p3 + DAYS_TO_1ST_AVTEMP * p4);
    tsq1 *= VarPar[30];
    DAYS_TO_1ST_SUMSTRS += p5 * (1.0 - WaterStress) + p6 * (1.0 - NStressVeg);
    tsq1 -= DAYS_TO_1ST_SUMSTRS;
    tsq1
}

unsafe fn create_first_square(stem_n_ratio: f64) {
    FruitingCode[0][0][0] = 1;
    FruitFraction[0][0][0] = 1.0;

    LeafAreaNodes[0][0][0] = VarPar[34];
    LeafWeightNodes[0][0][0] = VarPar[34] * LeafWeightAreaRatio;
    TotalStemWeight -= LeafWeightNodes[0][0][0];
    LeafNitrogen += LeafWeightNodes[0][0][0] * stem_n_ratio;
    StemNitrogen -= LeafWeightNodes[0][0][0] * stem_n_ratio;

    NumFruitBranches[0] = 1;
    NumNodes[0][0] = 1;
    NodeLayer[0][0] = ((PlantHeight / 5.0) as i32).clamp(0, 19);
    FruitGrowthRatio = 1.0;
    AvrgNodeTemper[0][0][0] = AvrgDailyTemp;

    let cotylwt = 0.20;
    AbscisedLeafWeight += cotylwt;
    CumPlantNLoss += cotylwt * LeafNitrogen / TotalLeafWeight();
    LeafNitrogen -= cotylwt * LeafNitrogen / TotalLeafWeight();
    PixInPlants -= cotylwt * pixcon;
}

unsafe fn add_vegetative_branch(
    delay_veg_by_c_stress: f64,
    stem_n_ratio: f64,
    days_to_1st_square: f64,
) {
    const VPVEGB: [f64; 3] = [13.39, -0.696, 0.012];

    let prev_branch = (NumVegBranches - 1) as usize;
    let time_to_next_veg_branch = VPVEGB[0]
        + AvrgNodeTemper[prev_branch][0][0]
            * (VPVEGB[1] + AvrgNodeTemper[prev_branch][0][0] * VPVEGB[2]);

    let threshold = time_to_next_veg_branch
        + delay_veg_by_c_stress
        + PHEN_DELAY_BY_N_STRESS
        + days_to_1st_square;
    if AgeOfSite[prev_branch][0][0] < threshold {
        return;
    }

    NumVegBranches += 1;
    if NumVegBranches > 3 {
        NumVegBranches = 3;
        return;
    }

    let new_branch = (NumVegBranches - 1) as usize;
    FruitFraction[new_branch][0][0] = 1.0;
    FruitingCode[new_branch][0][0] = 1;

    LeafAreaNodes[new_branch][0][0] = VarPar[34];
    LeafWeightNodes[new_branch][0][0] = VarPar[34] * LeafWeightAreaRatio;
    LeafAreaMainStem[new_branch][0] = VarPar[34];
    LeafWeightMainStem[new_branch][0] = LeafAreaMainStem[new_branch][0] * LeafWeightAreaRatio;

    TotalStemWeight -= LeafWeightNodes[new_branch][0][0] + LeafWeightMainStem[new_branch][0];
    let addlfn =
        (LeafWeightNodes[new_branch][0][0] + LeafWeightMainStem[new_branch][0]) * stem_n_ratio;
    LeafNitrogen += addlfn;
    StemNitrogen -= addlfn;

    AvrgNodeTemper[new_branch][0][0] = AvrgDailyTemp;
    NumFruitBranches[new_branch] = 1;
    NumNodes[new_branch][0] = 1;
    NodeLayer[new_branch][0] = ((PlantHeight / 5.0) as i32).clamp(0, 19);
}

unsafe fn add_fruiting_branch(k: usize, delay_veg_by_c_stress: f64, stem_n_ratio: f64) {
    const VFRTBR: [f64; 8] = [0.8, 0.95, 33.0, 4.461, -0.1912, 0.00265, 1.8, -1.32];

    DelayNewFruBranch[k] += delay_veg_by_c_stress + VFRTBR[0] * PHEN_DELAY_BY_N_STRESS;
    DelayNewFruBranch[k] += VFRTBR[1] * (1.0 - WaterStress);

    let nbrch = (NumFruitBranches[k] - 1) as usize;
    let mut tav = AvrgNodeTemper[k][nbrch][0];
    if tav > VFRTBR[2] {
        tav = VFRTBR[2];
    }

    let mut time_to_next_fru_branch =
        VarPar[35] + tav * (VFRTBR[3] + tav * (VFRTBR[4] + tav * VFRTBR[5]));
    if k > 0 {
        time_to_next_fru_branch *= VFRTBR[6];
    }
    time_to_next_fru_branch =
        time_to_next_fru_branch * (1.0 + VFRTBR[7] * (1.0 - DensityFactor)) + DelayNewFruBranch[k];

    if Kday > KdayAdjust && Kday <= (KdayAdjust + NumAdjustDays) && nadj[0] {
        if AdjAddMSNodesRate == 0.0 {
            time_to_next_fru_branch = 100.0;
        } else {
            time_to_next_fru_branch /= AdjAddMSNodesRate;
        }
    }

    if AgeOfSite[k][nbrch][0] < time_to_next_fru_branch {
        return;
    }

    NumFruitBranches[k] += 1;
    if NumFruitBranches[k] > 30 {
        NumFruitBranches[k] = 30;
        return;
    }

    let newbr = (NumFruitBranches[k] - 1) as usize;
    NumNodes[k][newbr] = 1;
    FruitFraction[k][newbr][0] = 1.0;
    FruitingCode[k][newbr][0] = 1;

    LeafAreaNodes[k][newbr][0] = VarPar[34];
    LeafWeightNodes[k][newbr][0] = VarPar[34] * LeafWeightAreaRatio;
    LeafAreaMainStem[k][newbr] = VarPar[34];
    LeafWeightMainStem[k][newbr] = LeafAreaMainStem[k][newbr] * LeafWeightAreaRatio;

    TotalStemWeight -= LeafWeightMainStem[k][newbr] + LeafWeightNodes[k][newbr][0];
    let addlfn = (LeafWeightMainStem[k][newbr] + LeafWeightNodes[k][newbr][0]) * stem_n_ratio;
    LeafNitrogen += addlfn;
    NodeLayer[k][newbr] = ((PlantHeight / 5.0) as i32).clamp(0, 19);
    StemNitrogen -= addlfn;

    AvrgNodeTemper[k][newbr][0] = AvrgDailyTemp;
    DelayNewFruBranch[k] = 0.0;
}

unsafe fn add_fruiting_node(k: usize, l: usize, delay_frt_by_c_stress: f64, stem_n_ratio: f64) {
    const VFRTNOD: [f64; 6] = [1.32, 0.90, 33.0, 7.6725, -0.3297, 0.004657];

    DelayNewNode[k][l] += (delay_frt_by_c_stress + VFRTNOD[0] * PHEN_DELAY_BY_N_STRESS) / pixdn;
    DelayNewNode[k][l] += VFRTNOD[1] * (1.0 - WaterStress);

    let nnid = (NumNodes[k][l] - 1) as usize;
    let mut tav = AvrgNodeTemper[k][l][nnid];
    if tav > VFRTNOD[2] {
        tav = VFRTNOD[2];
    }

    let mut time_to_next_fru_node =
        VarPar[36] + tav * (VFRTNOD[3] + tav * (VFRTNOD[4] + tav * VFRTNOD[5]));
    time_to_next_fru_node =
        time_to_next_fru_node * (1.0 + VarPar[37] * (1.0 - DensityFactor)) + DelayNewNode[k][l];

    if Kday > KdayAdjust && Kday <= (KdayAdjust + NumAdjustDays) && nadj[2] {
        if AdjAddSitesRate == 0.0 {
            time_to_next_fru_node = 100.0;
        } else {
            time_to_next_fru_node /= AdjAddSitesRate;
        }
    }

    if AgeOfSite[k][l][nnid] < time_to_next_fru_node {
        return;
    }

    NumNodes[k][l] += 1;
    if NumNodes[k][l] > 5 {
        NumNodes[k][l] = 5;
        return;
    }

    let newnod = nnid + 1;
    FruitFraction[k][l][newnod] = 1.0;
    FruitingCode[k][l][newnod] = 1;
    LeafAreaNodes[k][l][newnod] = VarPar[34];
    LeafWeightNodes[k][l][newnod] = VarPar[34] * LeafWeightAreaRatio;
    TotalStemWeight -= LeafWeightNodes[k][l][newnod];
    LeafNitrogen += LeafWeightNodes[k][l][newnod] * stem_n_ratio;
    StemNitrogen -= LeafWeightNodes[k][l][newnod] * stem_n_ratio;

    AvrgNodeTemper[k][l][newnod] = AvrgDailyTemp;
    DelayNewNode[k][l] = 0.0;
}

unsafe fn fruiting_site(k: usize, l: usize, m: usize, node_recent_white_flower: &mut i32) {
    const VFRSITE: [f64; 15] = [
        0.60, 0.40, 12.25, 0.40, 33.0, 0.20, 0.04, 0.45, 26.10, 9.0, 0.10, 3.0, 1.129, 0.043, 0.26,
    ];

    if FruitingCode[k][l][m] <= 0 {
        BOLTMP[k][l][m] = 0.0;
        return;
    }

    if NumFruitSites <= 0 {
        FIB_LENGTH = 0.0;
        FIB_STRENGTH = 0.0;
    }
    NumFruitSites += 1;

    let agefac = (1.0 - WaterStress) * VFRSITE[0] + (1.0 - NStressVeg) * VFRSITE[1];
    LeafAge[k][l][m] += DayInc + agefac;
    if DayFirstDef > 0 && Daynum > DayFirstDef {
        LeafAge[k][l][m] += VarPar[38];
    }

    if FruitingCode[k][l][m] >= 3 && FruitingCode[k][l][m] <= 6 {
        return;
    }

    let tmin = GetFromClim(CLIMATE_METRIC_TMIN, Daynum);
    let tmax = GetFromClim(CLIMATE_METRIC_TMAX, Daynum);
    let mut ageinc = DayInc;
    if tmin < VFRSITE[2] {
        ageinc += VFRSITE[3] * (VFRSITE[2] - tmin);
    }
    if tmax > VFRSITE[4] {
        let mut adjust = VFRSITE[6] * (tmax - VFRSITE[4]);
        if adjust > VFRSITE[5] {
            adjust = VFRSITE[5];
        }
        ageinc -= adjust;
    }

    if ageinc < VFRSITE[7] {
        ageinc = VFRSITE[7];
    }

    AvrgNodeTemper[k][l][m] = (AvrgNodeTemper[k][l][m] * AgeOfSite[k][l][m]
        + AvrgDailyTemp * ageinc)
        / (AgeOfSite[k][l][m] + ageinc);
    AgeOfSite[k][l][m] += ageinc;

    if FruitingCode[k][l][m] == 1 {
        if AgeOfSite[k][l][m] >= VFRSITE[8] {
            BOLTMP[k][l][m] = AvrgDailyTemp;
            AgeOfBoll[k][l][m] = DayInc;
            FruitingCode[k][l][m] = 7;
            new_boll_formation(k, l, m);
            if CottonWeightGreenBolls > 0.0 && FirstBloom <= 1 {
                FirstBloom = Daynum;
            }
            if k == 0 && m == 0 {
                if l as i32 > *node_recent_white_flower {
                    *node_recent_white_flower = l as i32;
                }
            }
        }
        return;
    }

    if BollWeight[k][l][m] > 0.0 {
        let dum = if LeafAreaIndex <= VFRSITE[11] && Kday > 100 {
            VFRSITE[12] - VFRSITE[13] * LeafAreaIndex
        } else {
            1.0
        };
        let dagebol = DayInc * dum
            + VFRSITE[14] * (1.0 - WaterStress)
            + VFRSITE[10] * (1.0 - NStressFruiting);
        BOLTMP[k][l][m] = (BOLTMP[k][l][m] * AgeOfBoll[k][l][m] + AvrgDailyTemp * dagebol)
            / (AgeOfBoll[k][l][m] + dagebol);
        AgeOfBoll[k][l][m] += dagebol;
    }

    if FruitingCode[k][l][m] == 7 {
        if AgeOfBoll[k][l][m] >= VFRSITE[9] {
            FruitingCode[k][l][m] = 2;
        }
        return;
    }

    if FruitingCode[k][l][m] == 2 {
        boll_opening(k, l, m, BOLTMP[k][l][m]);
    }
}

unsafe fn new_boll_formation(k: usize, l: usize, m: usize) {
    const SEED_RATIO: f64 = 0.64;
    const VNEWBOLL: [f64; 2] = [0.31, 0.02];

    if !bPollinSwitch {
        FruitingCode[k][l][m] = 6;
        FruitFraction[k][l][m] = 0.0;
        BloomWeightLoss += SquareWeight[k][l][m];
        SquareWeight[k][l][m] = 0.0;
        return;
    }

    let bolinit = VNEWBOLL[0] * SquareWeight[k][l][m];
    BollWeight[k][l][m] = 0.2 * bolinit;
    BurrWeight[k][l][m] = bolinit - BollWeight[k][l][m];
    BloomWeightLoss += SquareWeight[k][l][m] - bolinit;

    let mut sqr1n = SquareNConc * SquareWeight[k][l][m];
    SquareNitrogen -= sqr1n;
    CumPlantNLoss += sqr1n * (1.0 - VNEWBOLL[0]);
    sqr1n *= VNEWBOLL[0];

    let mut seed1n = BollWeight[k][l][m] * SEED_RATIO * VNEWBOLL[1];
    if seed1n > sqr1n {
        seed1n = sqr1n;
    }
    SeedNitrogen += seed1n;
    BurrNitrogen += sqr1n - seed1n;

    CottonWeightGreenBolls += BollWeight[k][l][m];
    BurrWeightGreenBolls += BurrWeight[k][l][m];
    TotalSquareWeight -= SquareWeight[k][l][m];
    SquareWeight[k][l][m] = 0.0;
}

unsafe fn boll_opening(k: usize, l: usize, m: usize, tmpboll: f64) {
    let ddpar1 = 1.0;
    let ddpar2 = 0.8;
    let vboldhs: [f64; 11] = [
        30.0, 41.189, -1.6057, 0.020743, 70.0, 0.994, 56.603, -2.921, 0.059, 1.219, 0.0065,
    ];

    let mut atn = tmpboll;
    if atn > vboldhs[0] {
        atn = vboldhs[0];
    }

    let mut dehiss = VarPar[39] + atn * (vboldhs[1] + atn * (vboldhs[2] + atn * vboldhs[3]));
    dehiss *= VarPar[40];
    if dehiss > vboldhs[4] {
        dehiss = vboldhs[4];
    }

    if DayFirstDef > 0 && Daynum > DayFirstDef {
        dehiss *= vboldhs[5].powf((Daynum - DayFirstDef) as f64);
    }

    if LeafAreaIndex < ddpar1 {
        let mut fdhslai = ddpar2 + LeafAreaIndex * (1.0 - ddpar2) / ddpar1;
        fdhslai = fdhslai.clamp(0.0, 1.0);
        dehiss *= fdhslai;
    }

    if AgeOfBoll[k][l][m] < dehiss {
        return;
    }

    FruitingCode[k][l][m] = 3;
    CottonWeightOpenBolls += BollWeight[k][l][m];
    BurrWeightOpenBolls += BurrWeight[k][l][m];
    CottonWeightGreenBolls -= BollWeight[k][l][m];
    BurrWeightGreenBolls -= BurrWeight[k][l][m];

    ginp = (VarPar[41] - VarPar[42] * atn) / 100.0;
    Gintot = (Gintot * NumOpenBolls + ginp * FruitFraction[k][l][m])
        / (NumOpenBolls + FruitFraction[k][l][m]);
    LintYield += ginp * BollWeight[k][l][m] * PlantPopulation * 0.001;

    let fsx = vboldhs[6] + atn * (vboldhs[7] + vboldhs[8] * atn);
    let flx = vboldhs[9] - vboldhs[10] * atn;
    FIB_STRENGTH = (FIB_STRENGTH * NumOpenBolls + fsx * FruitFraction[k][l][m])
        / (NumOpenBolls + FruitFraction[k][l][m]);
    FIB_LENGTH = (FIB_LENGTH * NumOpenBolls + flx * FruitFraction[k][l][m])
        / (NumOpenBolls + FruitFraction[k][l][m]);

    NumOpenBolls += FruitFraction[k][l][m];
}
