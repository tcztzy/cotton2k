use crate::general_functions::GetFromClim;
use crate::LegacyGlobalState;
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
    let mut legacy = LegacyGlobalState::from_globals();
    if legacy.leaf_area_index <= 0.0001 {
        return;
    }

    let droplf = drop_leaf_age(legacy.leaf_area_index);
    legacy.write_to_globals();
    pre_fruit_leaf_abscission(droplf);

    legacy = LegacyGlobalState::from_globals();
    for k in 0..legacy.num_veg_branches as usize {
        let nbrch = legacy.num_fruit_branches[k];
        for l in 0..nbrch as usize {
            main_stem_leaf_abscission(k, l, droplf);
        }
        legacy = LegacyGlobalState::from_globals();
    }

    if legacy.day_first_def > 0 && legacy.daynum >= legacy.day_first_def {
        legacy.write_to_globals();
        defoliation_leaf_abscission();
        legacy = LegacyGlobalState::from_globals();
    }

    if legacy.reserve_c > 0.0 {
        let resratio = 0.20;
        let resmax = resratio * TotalLeafWeight();
        if legacy.reserve_c > resmax {
            legacy.abscised_leaf_weight += legacy.reserve_c - resmax;
            legacy.reserve_c = resmax;
        }
    }

    legacy.leaf_area_index = TotalLeafArea() / legacy.per_plant_area;
    if legacy.leaf_area_index < 0.0001 {
        legacy.leaf_area_index = 0.0001;
    }
    legacy.write_to_globals();
}

fn pre_fruit_leaf_abscission(droplf: f64) {
    let mut legacy = LegacyGlobalState::from_globals();
    for j in 0..legacy.num_pre_fru_nodes as usize {
        if legacy.first_square > 0 {
            legacy.age_of_pre_fru_node[j] += legacy.day_inc;
        }
        if legacy.age_of_pre_fru_node[j] >= droplf
            && legacy.leaf_area_pre_fru[j] > 0.0
            && legacy.leaf_area_index > 0.1
        {
            legacy.leaf_area[legacy.node_layer_pre_fru[j] as usize] -= legacy.leaf_area_pre_fru[j];
            legacy.abscised_leaf_weight +=
                legacy.leaf_weight_pre_fru[j] + legacy.petiole_weight_pre_fru[j];
            legacy.total_petiole_weight -= legacy.petiole_weight_pre_fru[j];
            legacy.pix_in_plants -= legacy.leaf_weight_pre_fru[j] * legacy.pixcon;
            legacy.leaf_nitrogen -= legacy.leaf_weight_pre_fru[j] * legacy.leaf_n_conc;
            legacy.petiole_nitrogen -= legacy.petiole_weight_pre_fru[j] * legacy.petiole_n_conc;
            legacy.cum_plant_n_loss += legacy.leaf_weight_pre_fru[j] * legacy.leaf_n_conc
                + legacy.petiole_weight_pre_fru[j] * legacy.petiole_n_conc;
            legacy.leaf_area_pre_fru[j] = 0.0;
            legacy.leaf_weight_pre_fru[j] = 0.0;
            legacy.petiole_weight_pre_fru[j] = 0.0;
            if legacy.day_first_def > 0 && legacy.daynum > legacy.day_first_def {
                legacy.num_abscised_leaves += 1;
            }
        }
    }
    legacy.write_to_globals();
}

fn main_stem_leaf_abscission(k: usize, l: usize, droplf: f64) {
    let mut legacy = LegacyGlobalState::from_globals();
    if legacy.leaf_age[[k, l, 0]] > droplf
        && legacy.leaf_area_main_stem[[k, l]] > 0.0
        && legacy.leaf_area_index > 0.1
    {
        legacy.abscised_leaf_weight +=
            legacy.leaf_weight_main_stem[[k, l]] + legacy.petiole_weight_main_stem[[k, l]];
        legacy.total_petiole_weight -= legacy.petiole_weight_main_stem[[k, l]];
        legacy.pix_in_plants -= legacy.leaf_weight_main_stem[[k, l]] * legacy.pixcon;
        legacy.leaf_nitrogen -= legacy.leaf_weight_main_stem[[k, l]] * legacy.leaf_n_conc;
        legacy.petiole_nitrogen -= legacy.petiole_weight_main_stem[[k, l]] * legacy.petiole_n_conc;
        legacy.cum_plant_n_loss += legacy.leaf_weight_main_stem[[k, l]] * legacy.leaf_n_conc
            + legacy.petiole_weight_main_stem[[k, l]] * legacy.petiole_n_conc;
        legacy.leaf_area[legacy.node_layer[[k, l]] as usize] -= legacy.leaf_area_main_stem[[k, l]];
        legacy.leaf_area_main_stem[[k, l]] = 0.0;
        legacy.leaf_weight_main_stem[[k, l]] = 0.0;
        legacy.petiole_weight_main_stem[[k, l]] = 0.0;
        if legacy.day_first_def > 0 && legacy.daynum > legacy.day_first_def {
            legacy.num_abscised_leaves += 1;
        }
    }

    let nnid = legacy.num_nodes[[k, l]];
    legacy.write_to_globals();
    for m in 0..nnid as usize {
        fruit_node_leaf_abscission(k, l, m, droplf);
    }
}

fn fruit_node_leaf_abscission(k: usize, l: usize, m: usize, droplf: f64) {
    let mut legacy = LegacyGlobalState::from_globals();
    if legacy.leaf_age[[k, l, m]] >= droplf
        && legacy.leaf_area_nodes[[k, l, m]] > 0.0
        && legacy.leaf_area_index > 0.1
    {
        legacy.abscised_leaf_weight +=
            legacy.leaf_weight_nodes[[k, l, m]] + legacy.petiole_weight_nodes[[k, l, m]];
        legacy.total_petiole_weight -= legacy.petiole_weight_nodes[[k, l, m]];
        legacy.pix_in_plants -= legacy.leaf_weight_nodes[[k, l, m]] * legacy.pixcon;
        legacy.leaf_nitrogen -= legacy.leaf_weight_nodes[[k, l, m]] * legacy.leaf_n_conc;
        legacy.petiole_nitrogen -= legacy.petiole_weight_nodes[[k, l, m]] * legacy.petiole_n_conc;
        legacy.cum_plant_n_loss += legacy.leaf_weight_nodes[[k, l, m]] * legacy.leaf_n_conc
            + legacy.petiole_weight_nodes[[k, l, m]] * legacy.petiole_n_conc;
        legacy.leaf_area[legacy.node_layer[[k, l]] as usize] -= legacy.leaf_area_nodes[[k, l, m]];
        legacy.leaf_area_nodes[[k, l, m]] = 0.0;
        legacy.leaf_weight_nodes[[k, l, m]] = 0.0;
        legacy.petiole_weight_nodes[[k, l, m]] = 0.0;
        if legacy.day_first_def > 0 && legacy.daynum > legacy.day_first_def {
            legacy.num_abscised_leaves += 1;
        }
    }
    legacy.write_to_globals();
}

fn defoliation_leaf_abscission() {
    let mut legacy = LegacyGlobalState::from_globals();
    if legacy.daynum == legacy.day_first_def {
        for j in 0..legacy.num_pre_fru_nodes as usize {
            if legacy.leaf_area_pre_fru[j] > 0.0 {
                legacy.leaf_area[NodeLayerPreFru(j)] -= legacy.leaf_area_pre_fru[j];
                legacy.leaf_area_pre_fru[j] = 0.0;
                legacy.abscised_leaf_weight +=
                    legacy.leaf_weight_pre_fru[j] + legacy.petiole_weight_pre_fru[j];
                legacy.total_petiole_weight -= legacy.petiole_weight_pre_fru[j];
                legacy.pix_in_plants -= legacy.leaf_weight_pre_fru[j] * legacy.pixcon;
                legacy.leaf_nitrogen -= legacy.leaf_weight_pre_fru[j] * legacy.leaf_n_conc;
                legacy.petiole_nitrogen -= legacy.petiole_weight_pre_fru[j] * legacy.petiole_n_conc;
                legacy.cum_plant_n_loss += legacy.leaf_weight_pre_fru[j] * legacy.leaf_n_conc
                    + legacy.petiole_weight_pre_fru[j] * legacy.petiole_n_conc;
                legacy.leaf_weight_pre_fru[j] = 0.0;
                legacy.petiole_weight_pre_fru[j] = 0.0;
            }
        }
    }

    if legacy.daynum <= legacy.day_first_def {
        legacy.write_to_globals();
        return;
    }

    let mut sort_by_age = [0.0; 450];
    let mut indexk = [0; 450];
    let mut indexl = [0; 450];
    let mut indexm = [0; 450];

    let mut lefcnt: usize = 0;
    for k in 0..legacy.num_veg_branches as usize {
        let nbrch = legacy.num_fruit_branches[k];
        for l in 0..nbrch as usize {
            if legacy.leaf_weight_main_stem[[k, l]] > 0.0 {
                sort_by_age[lefcnt] = legacy.age_of_site[[k, l, 0]];
                indexk[lefcnt] = k as i32;
                indexl[lefcnt] = l as i32;
                indexm[lefcnt] = 66;
                lefcnt += 1;
            }
            let nnid = legacy.num_nodes[[k, l]];
            for m in 0..nnid as usize {
                if legacy.leaf_weight_nodes[[k, l, m]] > 0.0 {
                    sort_by_age[lefcnt] = legacy.age_of_site[[k, l, m]];
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

    let mut num_leaves_to_shed = (lefcnt as f64 * legacy.percent_defoliation / 100.0) as i32;
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
                legacy.abscised_leaf_weight +=
                    legacy.leaf_weight_main_stem[[k, l]] + legacy.petiole_weight_main_stem[[k, l]];
                legacy.total_petiole_weight -= legacy.petiole_weight_main_stem[[k, l]];
                pixlos = legacy.leaf_weight_main_stem[[k, l]] * legacy.pixcon;
                legacy.leaf_nitrogen -= legacy.leaf_weight_main_stem[[k, l]] * legacy.leaf_n_conc;
                legacy.petiole_nitrogen -=
                    legacy.petiole_weight_main_stem[[k, l]] * legacy.petiole_n_conc;
                legacy.cum_plant_n_loss += legacy.leaf_weight_main_stem[[k, l]]
                    * legacy.leaf_n_conc
                    + legacy.petiole_weight_main_stem[[k, l]] * legacy.petiole_n_conc;
                legacy.leaf_area[legacy.node_layer[[k, l]] as usize] -=
                    legacy.leaf_area_main_stem[[k, l]];
                legacy.leaf_area_main_stem[[k, l]] = 0.0;
                legacy.leaf_weight_main_stem[[k, l]] = 0.0;
                legacy.petiole_weight_main_stem[[k, l]] = 0.0;
            } else {
                let m = m as usize;
                legacy.abscised_leaf_weight +=
                    legacy.leaf_weight_nodes[[k, l, m]] + legacy.petiole_weight_nodes[[k, l, m]];
                legacy.total_petiole_weight -= legacy.petiole_weight_nodes[[k, l, m]];
                pixlos = legacy.leaf_weight_nodes[[k, l, m]] * legacy.pixcon;
                legacy.leaf_nitrogen -= legacy.leaf_weight_nodes[[k, l, m]] * legacy.leaf_n_conc;
                legacy.petiole_nitrogen -=
                    legacy.petiole_weight_nodes[[k, l, m]] * legacy.petiole_n_conc;
                legacy.cum_plant_n_loss += legacy.leaf_weight_nodes[[k, l, m]] * legacy.leaf_n_conc
                    + legacy.petiole_weight_nodes[[k, l, m]] * legacy.petiole_n_conc;
                legacy.leaf_area[legacy.node_layer[[k, l]] as usize] -=
                    legacy.leaf_area_nodes[[k, l, m]];
                legacy.leaf_area_nodes[[k, l, m]] = 0.0;
                legacy.leaf_weight_nodes[[k, l, m]] = 0.0;
                legacy.petiole_weight_nodes[[k, l, m]] = 0.0;
            }
            legacy.pix_in_plants -= pixlos;
            num_leaves_to_shed -= 1;
            legacy.num_abscised_leaves += 1;
        }
    }
    legacy.write_to_globals();
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
    const VABSFR: [f64; 9] = [21.0, 0.42, 30.0, 0.05, 6.0, 2.25, 0.60, 5.0, 0.20];

    let mut legacy = LegacyGlobalState::from_globals();
    legacy.num_shedding_tags += 1;
    if legacy.num_shedding_tags > 1 {
        for lt in (1..legacy.num_shedding_tags as usize).rev() {
            let ltm1 = lt - 1;
            legacy.shed_by_carbon_stress[lt] = legacy.shed_by_carbon_stress[ltm1];
            legacy.shed_by_nitrogen_stress[lt] = legacy.shed_by_nitrogen_stress[ltm1];
            legacy.shed_by_water_stress[lt] = legacy.shed_by_water_stress[ltm1];
            legacy.abscission_lag[lt] = legacy.abscission_lag[ltm1];
        }
    }

    if legacy.carbon_stress < legacy.var_par[43] {
        legacy.shed_by_carbon_stress[0] =
            (legacy.var_par[43] - legacy.carbon_stress) / legacy.var_par[43];
    } else {
        legacy.shed_by_carbon_stress[0] = 0.0;
    }
    if legacy.nitrogen_stress < VABSFR[1] {
        legacy.shed_by_nitrogen_stress[0] = (VABSFR[1] - legacy.nitrogen_stress) / VABSFR[1];
    } else {
        legacy.shed_by_nitrogen_stress[0] = 0.0;
    }
    if legacy.water_stress < legacy.var_par[44] {
        legacy.shed_by_water_stress[0] =
            (legacy.var_par[44] - legacy.water_stress) / legacy.var_par[44];
    } else {
        legacy.shed_by_water_stress[0] = 0.0;
    }
    legacy.abscission_lag[0] = 0.01;

    let tmax = GetFromClim(CLIMATE_METRIC_TMAX, legacy.daynum);
    for lt in 0..legacy.num_shedding_tags as usize {
        legacy.abscission_lag[lt] += legacy.day_inc.max(0.40);
        if tmax > VABSFR[2] {
            legacy.abscission_lag[lt] += (tmax - VABSFR[2]) * VABSFR[3];
        }
    }
    legacy.write_to_globals();

    let mut idecr = 0;
    let mut legacy = LegacyGlobalState::from_globals();
    for lt in 0..legacy.num_shedding_tags as usize {
        if legacy.abscission_lag[lt] >= VABSFR[4] || lt >= 20 {
            let gin1 = if legacy.gintot > 0.0 {
                legacy.gintot
            } else {
                legacy.ginp
            };
            for k in 0..legacy.num_veg_branches as usize {
                let nbrch = legacy.num_fruit_branches[k];
                for l in 0..nbrch as usize {
                    let nnid = legacy.num_nodes[[k, l]];
                    for m in 0..nnid as usize {
                        let code = legacy.fruiting_code[[k, l, m]];
                        if code == 1 || code == 2 || code == 7 {
                            let abscission_ratio = site_abscission_ratio(k, l, m, lt);
                            if abscission_ratio > 0.0 {
                                if code == 1 {
                                    square_abscission(k, l, m, abscission_ratio);
                                } else {
                                    boll_abscission(k, l, m, abscission_ratio, gin1);
                                }
                            }
                        }
                    }
                }
            }
            legacy = LegacyGlobalState::from_globals();
            legacy.shed_by_carbon_stress[lt] = 0.0;
            legacy.shed_by_nitrogen_stress[lt] = 0.0;
            legacy.shed_by_water_stress[lt] = 0.0;
            legacy.abscission_lag[lt] = 0.0;
            idecr += 1;
        }
    }

    legacy.num_shedding_tags -= idecr;
    let should_adjust = legacy.kday > legacy.kday_adjust
        && legacy.kday <= (legacy.kday_adjust + legacy.num_adjust_days);
    legacy.write_to_globals();
    if should_adjust {
        adjust_abscission();
    }

    compute_site_numbers();
}

fn site_abscission_ratio(k: usize, l: usize, m: usize, lt: usize) -> f64 {
    const VABSC: [f64; 5] = [21.0, 2.25, 0.60, 5.0, 0.20];

    let legacy = LegacyGlobalState::from_globals();
    let mut pabs = 0.0;
    let mut shedt = 0.0;
    if legacy.fruiting_code[[k, l, m]] == 1 {
        if legacy.age_of_site[[k, l, m]] < VABSC[3] {
            pabs = 0.0;
        } else {
            let xsqage = legacy.age_of_site[[k, l, m]] - VABSC[3];
            if xsqage >= VABSC[0] {
                pabs = legacy.var_par[46];
            } else {
                pabs = legacy.var_par[46]
                    + (legacy.var_par[45] - legacy.var_par[46])
                        * ((VABSC[0] - xsqage) / VABSC[0]).powf(VABSC[1]);
            }
        }
        shedt = 1.0
            - (1.0 - legacy.shed_by_carbon_stress[lt]) * (1.0 - legacy.shed_by_nitrogen_stress[lt]);
    } else if legacy.fruiting_code[[k, l, m]] == 7
        && legacy.age_of_boll[[k, l, m]] <= legacy.var_par[47]
    {
        pabs = legacy.var_par[48];
        shedt = 1.0
            - (1.0 - legacy.shed_by_carbon_stress[lt])
                * (1.0 - VABSC[2] * legacy.shed_by_nitrogen_stress[lt]);
    } else if legacy.age_of_boll[[k, l, m]] > legacy.var_par[47]
        && legacy.age_of_boll[[k, l, m]] <= (legacy.var_par[47] + legacy.var_par[49])
    {
        pabs = legacy.var_par[48]
            - (legacy.var_par[48] - legacy.var_par[50])
                * (legacy.age_of_boll[[k, l, m]] - legacy.var_par[47])
                / legacy.var_par[49];
        shedt = 1.0
            - (1.0 - legacy.shed_by_carbon_stress[lt])
                * (1.0 - VABSC[4] * legacy.shed_by_nitrogen_stress[lt])
                * (1.0 - legacy.shed_by_water_stress[lt]);
    } else if legacy.age_of_boll[[k, l, m]] > (legacy.var_par[47] + legacy.var_par[49])
        && legacy.age_of_boll[[k, l, m]] <= (legacy.var_par[47] + 2.0 * legacy.var_par[49])
    {
        pabs = legacy.var_par[50] / legacy.var_par[49]
            * (legacy.var_par[47] + 2.0 * legacy.var_par[49] - legacy.age_of_boll[[k, l, m]]);
        shedt = legacy.shed_by_water_stress[lt];
    } else if legacy.age_of_boll[[k, l, m]] > (legacy.var_par[47] + 2.0 * legacy.var_par[49]) {
        pabs = 0.0;
    }

    let mut abscission_ratio = pabs * shedt * legacy.day_inc;
    if abscission_ratio > 1.0 {
        abscission_ratio = 1.0;
    }
    abscission_ratio
}

fn square_abscission(k: usize, l: usize, m: usize, abscission_ratio: f64) {
    let mut legacy = LegacyGlobalState::from_globals();
    let wtlos = legacy.square_weight[[k, l, m]] * abscission_ratio;
    legacy.square_nitrogen -= wtlos * legacy.square_n_conc;
    legacy.cum_plant_n_loss += wtlos * legacy.square_n_conc;
    legacy.square_weight[[k, l, m]] -= wtlos;
    legacy.bloom_weight_loss += wtlos;
    legacy.total_square_weight -= wtlos;
    legacy.fruit_fraction[[k, l, m]] *= 1.0 - abscission_ratio;

    if legacy.fruit_fraction[[k, l, m]] <= 0.001 {
        legacy.fruit_fraction[[k, l, m]] = 0.0;
        legacy.square_nitrogen -= legacy.square_weight[[k, l, m]] * legacy.square_n_conc;
        legacy.cum_plant_n_loss += legacy.square_weight[[k, l, m]] * legacy.square_n_conc;
        legacy.bloom_weight_loss += legacy.square_weight[[k, l, m]];
        legacy.total_square_weight -= legacy.square_weight[[k, l, m]];
        legacy.square_weight[[k, l, m]] = 0.0;
        legacy.fruiting_code[[k, l, m]] = 5;
    }
    legacy.write_to_globals();
}

fn boll_abscission(k: usize, l: usize, m: usize, abscission_ratio: f64, gin1: f64) {
    let mut legacy = LegacyGlobalState::from_globals();
    legacy.seed_nitrogen -=
        legacy.boll_weight[[k, l, m]] * abscission_ratio * (1.0 - gin1) * legacy.seed_n_conc;
    legacy.burr_nitrogen -= legacy.burr_weight[[k, l, m]] * abscission_ratio * legacy.burr_n_conc;
    legacy.cum_plant_n_loss +=
        legacy.boll_weight[[k, l, m]] * abscission_ratio * (1.0 - gin1) * legacy.seed_n_conc;
    legacy.cum_plant_n_loss +=
        legacy.burr_weight[[k, l, m]] * abscission_ratio * legacy.burr_n_conc;
    legacy.pix_in_plants -= (legacy.boll_weight[[k, l, m]] + legacy.burr_weight[[k, l, m]])
        * abscission_ratio
        * legacy.pixcon;
    legacy.green_bolls_lost +=
        (legacy.boll_weight[[k, l, m]] + legacy.burr_weight[[k, l, m]]) * abscission_ratio;
    legacy.cotton_weight_green_bolls -= legacy.boll_weight[[k, l, m]] * abscission_ratio;
    legacy.burr_weight_green_bolls -= legacy.burr_weight[[k, l, m]] * abscission_ratio;
    legacy.boll_weight[[k, l, m]] -= legacy.boll_weight[[k, l, m]] * abscission_ratio;
    legacy.burr_weight[[k, l, m]] -= legacy.burr_weight[[k, l, m]] * abscission_ratio;
    legacy.fruit_fraction[[k, l, m]] -= legacy.fruit_fraction[[k, l, m]] * abscission_ratio;

    if legacy.fruit_fraction[[k, l, m]] <= 0.001 {
        legacy.fruiting_code[[k, l, m]] = 4;
        legacy.seed_nitrogen -= legacy.boll_weight[[k, l, m]] * (1.0 - gin1) * legacy.seed_n_conc;
        legacy.burr_nitrogen -= legacy.burr_weight[[k, l, m]] * legacy.burr_n_conc;
        legacy.cum_plant_n_loss +=
            legacy.boll_weight[[k, l, m]] * (1.0 - gin1) * legacy.seed_n_conc;
        legacy.cum_plant_n_loss += legacy.burr_weight[[k, l, m]] * legacy.burr_n_conc;
        legacy.pix_in_plants -=
            (legacy.boll_weight[[k, l, m]] + legacy.burr_weight[[k, l, m]]) * legacy.pixcon;
        legacy.fruit_fraction[[k, l, m]] = 0.0;
        legacy.cotton_weight_green_bolls -= legacy.boll_weight[[k, l, m]];
        legacy.burr_weight_green_bolls -= legacy.burr_weight[[k, l, m]];
        legacy.green_bolls_lost += legacy.boll_weight[[k, l, m]] + legacy.burr_weight[[k, l, m]];
        legacy.boll_weight[[k, l, m]] = 0.0;
        legacy.burr_weight[[k, l, m]] = 0.0;
    }
    legacy.write_to_globals();
}

fn adjust_abscission() {
    let mut legacy = LegacyGlobalState::from_globals();
    let mut jx = [0_i32; 2];
    let mut abscsq = 0.0;
    if legacy.nadj[3] && legacy.adj_square_absc > 0.0 {
        jx[0] = 1;
        abscsq = legacy.adj_square_absc;
    }

    let mut abscgb = 0.0;
    if legacy.nadj[4] && legacy.adj_green_boll_absc > 0.0 {
        jx[1] = 1;
        let mut gbolnum = 0.0;
        for k in 0..legacy.num_veg_branches as usize {
            let nbrch = legacy.num_fruit_branches[k];
            for l in 0..nbrch as usize {
                let nnid = legacy.num_nodes[[k, l]];
                for m in 0..nnid as usize {
                    if legacy.fruiting_code[[k, l, m]] == 7 {
                        gbolnum += legacy.fruit_fraction[[k, l, m]];
                    }
                }
            }
        }
        if gbolnum > 0.0 {
            abscgb = legacy.adj_green_boll_absc * legacy.num_green_bolls / gbolnum;
        }
    }

    let mut gin1 = 0.0;
    if jx[1] == 1 {
        gin1 = if legacy.gintot > 0.0 {
            legacy.gintot
        } else {
            legacy.ginp
        };
    }

    legacy.write_to_globals();
    legacy = LegacyGlobalState::from_globals();
    for k in 0..legacy.num_veg_branches as usize {
        let nbrch = legacy.num_fruit_branches[k];
        for l in 0..nbrch as usize {
            let nnid = legacy.num_nodes[[k, l]];
            for m in 0..nnid as usize {
                if jx[0] == 1 && legacy.fruiting_code[[k, l, m]] == 1 {
                    adjust_square_abscission(k, l, m, abscsq);
                }
                if jx[1] == 1 && legacy.fruiting_code[[k, l, m]] == 7 {
                    adjust_young_boll_abscission(k, l, m, abscgb, gin1);
                }
            }
        }
        legacy = LegacyGlobalState::from_globals();
    }
}

fn adjust_square_abscission(k: usize, l: usize, m: usize, abscsq: f64) {
    let mut legacy = LegacyGlobalState::from_globals();
    let mut wtlos = legacy.square_weight[[k, l, m]] * abscsq;
    legacy.square_nitrogen -= wtlos * legacy.square_n_conc;
    legacy.cum_plant_n_loss += wtlos * legacy.square_n_conc;
    legacy.square_weight[[k, l, m]] -= wtlos;
    legacy.bloom_weight_loss += wtlos;
    legacy.total_square_weight -= wtlos;
    legacy.fruit_fraction[[k, l, m]] *= 1.0 - abscsq;

    if legacy.fruit_fraction[[k, l, m]] <= 0.001 {
        legacy.fruit_fraction[[k, l, m]] = 0.0;
        legacy.square_nitrogen -= legacy.square_weight[[k, l, m]] * legacy.square_n_conc;
        legacy.cum_plant_n_loss += legacy.square_weight[[k, l, m]] * legacy.square_n_conc;
        legacy.bloom_weight_loss += legacy.square_weight[[k, l, m]];
        legacy.total_square_weight -= legacy.square_weight[[k, l, m]];
        legacy.square_weight[[k, l, m]] = 0.0;
        legacy.fruiting_code[[k, l, m]] = 5;
    }

    if legacy.fruit_fraction[[k, l, m]] > 1.0 {
        wtlos = legacy.square_weight[[k, l, m]] * (1.0 - 1.0 / legacy.fruit_fraction[[k, l, m]]);
        legacy.fruit_fraction[[k, l, m]] = 1.0;
        legacy.square_nitrogen -= wtlos * legacy.square_n_conc;
        legacy.cum_plant_n_loss += wtlos * legacy.square_n_conc;
        legacy.square_weight[[k, l, m]] -= wtlos;
        legacy.bloom_weight_loss += wtlos;
        legacy.total_square_weight -= wtlos;
    }
    legacy.write_to_globals();
}

fn adjust_young_boll_abscission(k: usize, l: usize, m: usize, abscgb: f64, gin1: f64) {
    let mut legacy = LegacyGlobalState::from_globals();
    legacy.seed_nitrogen -=
        legacy.boll_weight[[k, l, m]] * abscgb * (1.0 - gin1) * legacy.seed_n_conc;
    legacy.burr_nitrogen -= legacy.burr_weight[[k, l, m]] * abscgb * legacy.burr_n_conc;
    legacy.cum_plant_n_loss +=
        legacy.boll_weight[[k, l, m]] * abscgb * (1.0 - gin1) * legacy.seed_n_conc;
    legacy.cum_plant_n_loss += legacy.burr_weight[[k, l, m]] * abscgb * legacy.burr_n_conc;
    let pixlos =
        (legacy.boll_weight[[k, l, m]] + legacy.burr_weight[[k, l, m]]) * abscgb * legacy.pixcon;
    legacy.pix_in_plants -= pixlos;
    legacy.green_bolls_lost +=
        (legacy.boll_weight[[k, l, m]] + legacy.burr_weight[[k, l, m]]) * abscgb;
    legacy.cotton_weight_green_bolls -= legacy.boll_weight[[k, l, m]] * abscgb;
    legacy.burr_weight_green_bolls -= legacy.burr_weight[[k, l, m]] * abscgb;
    legacy.boll_weight[[k, l, m]] *= 1.0 - abscgb;
    legacy.burr_weight[[k, l, m]] *= 1.0 - abscgb;
    legacy.fruit_fraction[[k, l, m]] *= 1.0 - abscgb;
    legacy.write_to_globals();

    adjust_boll_abscission(k, l, m, 2, gin1);
}

fn adjust_boll_abscission(k: usize, l: usize, m: usize, jx: i32, gin1: f64) {
    let mut legacy = LegacyGlobalState::from_globals();
    if legacy.fruit_fraction[[k, l, m]] <= 0.001 {
        legacy.seed_nitrogen -= legacy.boll_weight[[k, l, m]] * (1.0 - gin1) * legacy.seed_n_conc;
        legacy.burr_nitrogen -= legacy.burr_weight[[k, l, m]] * legacy.burr_n_conc;
        legacy.cum_plant_n_loss +=
            legacy.boll_weight[[k, l, m]] * (1.0 - gin1) * legacy.seed_n_conc;
        legacy.cum_plant_n_loss += legacy.burr_weight[[k, l, m]] * legacy.burr_n_conc;
        legacy.pix_in_plants -=
            (legacy.boll_weight[[k, l, m]] + legacy.burr_weight[[k, l, m]]) * legacy.pixcon;
        legacy.green_bolls_lost += legacy.boll_weight[[k, l, m]] + legacy.burr_weight[[k, l, m]];

        if jx == 2 {
            legacy.cotton_weight_green_bolls -= legacy.boll_weight[[k, l, m]];
            legacy.burr_weight_green_bolls -= legacy.burr_weight[[k, l, m]];
        } else if jx == 3 {
            legacy.cotton_weight_open_bolls -= legacy.boll_weight[[k, l, m]];
            legacy.burr_weight_open_bolls -= legacy.burr_weight[[k, l, m]];
        }

        legacy.fruit_fraction[[k, l, m]] = 0.0;
        legacy.boll_weight[[k, l, m]] = 0.0;
        legacy.burr_weight[[k, l, m]] = 0.0;
        legacy.fruiting_code[[k, l, m]] = 4;
    }

    if legacy.fruit_fraction[[k, l, m]] > 1.0 {
        let bolwlos =
            legacy.boll_weight[[k, l, m]] * (1.0 - 1.0 / legacy.fruit_fraction[[k, l, m]]);
        let burwlos =
            legacy.burr_weight[[k, l, m]] * (1.0 - 1.0 / legacy.fruit_fraction[[k, l, m]]);
        legacy.fruit_fraction[[k, l, m]] = 1.0;
        legacy.seed_nitrogen -= bolwlos * (1.0 - gin1) * legacy.seed_n_conc;
        legacy.burr_nitrogen -= burwlos * legacy.burr_n_conc;
        legacy.cum_plant_n_loss += bolwlos * (1.0 - gin1) * legacy.seed_n_conc;
        legacy.cum_plant_n_loss += burwlos * legacy.burr_n_conc;
        legacy.pix_in_plants -= (bolwlos + burwlos) * legacy.pixcon;
        legacy.green_bolls_lost += bolwlos + burwlos;
        legacy.boll_weight[[k, l, m]] -= bolwlos;
        legacy.burr_weight[[k, l, m]] -= burwlos;

        if jx == 2 {
            legacy.cotton_weight_green_bolls -= bolwlos;
            legacy.burr_weight_green_bolls -= burwlos;
        } else if jx == 3 {
            legacy.cotton_weight_open_bolls -= bolwlos;
            legacy.burr_weight_open_bolls -= burwlos;
        }
    }
    legacy.write_to_globals();
}

fn compute_site_numbers() {
    let mut legacy = LegacyGlobalState::from_globals();
    legacy.num_squares = 0.0;
    legacy.num_green_bolls = 0.0;
    legacy.num_open_bolls = 0.0;
    for k in 0..legacy.num_veg_branches as usize {
        let nbrch = legacy.num_fruit_branches[k];
        for l in 0..nbrch as usize {
            let nnid = legacy.num_nodes[[k, l]];
            for m in 0..nnid as usize {
                if legacy.fruiting_code[[k, l, m]] == 1 {
                    legacy.num_squares += legacy.fruit_fraction[[k, l, m]];
                } else if legacy.fruiting_code[[k, l, m]] == 7
                    || legacy.fruiting_code[[k, l, m]] == 2
                {
                    legacy.num_green_bolls += legacy.fruit_fraction[[k, l, m]];
                } else if legacy.fruiting_code[[k, l, m]] == 3 {
                    legacy.num_open_bolls += legacy.fruit_fraction[[k, l, m]];
                }
            }
        }
    }

    legacy.abscised_fruit_sites = legacy.num_fruit_sites as f64
        - legacy.num_squares
        - legacy.num_green_bolls
        - legacy.num_open_bolls;
    legacy.write_to_globals();
}

#[inline]
fn NodeLayerPreFru(index: usize) -> usize {
    LegacyGlobalState::from_globals().node_layer_pre_fru[index] as usize
}
