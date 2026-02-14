use crate::general_functions::GetFromClim;
use crate::plant::abscission::{fruiting_sites_abscission, leaf_abscission};
use crate::{TotalLeafWeight, CLIMATE_METRIC_TMAX, CLIMATE_METRIC_TMIN};
use std::sync::{LazyLock, RwLock};

#[derive(Debug, Clone, Copy)]
struct PhenologyScratch {
    fib_length: f64,
    fib_strength: f64,
    phen_delay_by_n_stress: f64,
    nwfl: i32,
    days_to_1st_square: f64,
    boltmp: [[[f64; 5]; 30]; 3],
    days_to_1st_avtemp: f64,
    days_to_1st_sumstrs: f64,
}

static PHENOLOGY_SCRATCH: LazyLock<RwLock<PhenologyScratch>> = LazyLock::new(|| {
    RwLock::new(PhenologyScratch {
        fib_length: 0.0,
        fib_strength: 0.0,
        phen_delay_by_n_stress: 0.0,
        nwfl: 0,
        days_to_1st_square: 0.0,
        boltmp: [[[0.0; 5]; 30]; 3],
        days_to_1st_avtemp: 0.0,
        days_to_1st_sumstrs: 0.0,
    })
});

fn read_phenology_scratch() -> PhenologyScratch {
    *PHENOLOGY_SCRATCH
        .read()
        .expect("phenology scratch state lock should not be poisoned")
}

fn with_phenology_scratch_mut<R>(f: impl FnOnce(&mut PhenologyScratch) -> R) -> R {
    let mut scratch = PHENOLOGY_SCRATCH
        .write()
        .expect("phenology scratch state lock should not be poisoned");
    f(&mut scratch)
}

pub fn cotton_phenology() {
    let mut scratch = read_phenology_scratch();
    const VPHENO: [f64; 8] = [0.65, -0.83, -1.67, -0.25, -0.75, 10.0, 15.0, 7.10];

    let mut legacy = crate::LegacyGlobalState::from_globals();
    legacy.num_fruit_sites = 0;
    let stem_n_ratio = legacy.stem_nitrogen / legacy.total_stem_weight;

    scratch.phen_delay_by_n_stress = VPHENO[0] * (1.0 - legacy.n_stress_veg);
    scratch.phen_delay_by_n_stress = scratch.phen_delay_by_n_stress.clamp(0.0, 1.0);

    let mut delay_veg_by_c_stress =
        legacy.var_par[27] + legacy.carbon_stress * (VPHENO[3] + VPHENO[4] * legacy.carbon_stress);
    delay_veg_by_c_stress = delay_veg_by_c_stress.clamp(0.0, 1.0);

    let mut delay_frt_by_c_stress =
        legacy.var_par[28] + legacy.carbon_stress * (VPHENO[1] + VPHENO[2] * legacy.carbon_stress);
    if delay_frt_by_c_stress > legacy.var_par[29] {
        delay_frt_by_c_stress = legacy.var_par[29];
    }
    if delay_frt_by_c_stress < 0.0 {
        delay_frt_by_c_stress = 0.0;
    }

    if legacy.first_square <= 0 {
        scratch.days_to_1st_square = days_to_first_square();
        legacy.write_to_globals();
        pre_fruiting_node(stem_n_ratio);
        legacy = crate::LegacyGlobalState::from_globals();
        if legacy.kday >= scratch.days_to_1st_square as i32 {
            legacy.first_square = legacy.daynum;
            legacy.write_to_globals();
            create_first_square(stem_n_ratio);
            legacy = crate::LegacyGlobalState::from_globals();
        } else {
            legacy.write_to_globals();
            leaf_abscission();
            with_phenology_scratch_mut(|state| *state = scratch);
            return;
        }
    } else {
        legacy.write_to_globals();
    }

    legacy = crate::LegacyGlobalState::from_globals();
    if legacy.num_veg_branches == 1 && legacy.per_plant_area >= VPHENO[5] {
        add_vegetative_branch(
            delay_veg_by_c_stress,
            stem_n_ratio,
            scratch.days_to_1st_square,
        );
    }
    legacy = crate::LegacyGlobalState::from_globals();
    if legacy.num_veg_branches == 2 && legacy.per_plant_area >= VPHENO[6] {
        add_vegetative_branch(
            delay_veg_by_c_stress,
            stem_n_ratio,
            scratch.days_to_1st_square,
        );
    }

    legacy = crate::LegacyGlobalState::from_globals();
    let mut nidmax = (VPHENO[7] * legacy.density_factor + 0.5) as i32;
    if nidmax > 5 {
        nidmax = 5;
    }

    let mut nwfl = scratch.nwfl;
    for k in 0..legacy.num_veg_branches as usize {
        if legacy.num_fruit_branches[k] < 30 {
            add_fruiting_branch(k, delay_veg_by_c_stress, stem_n_ratio);
        }
        legacy = crate::LegacyGlobalState::from_globals();
        for l in 0..legacy.num_fruit_branches[k] as usize {
            if legacy.num_nodes[[k, l]] < nidmax {
                add_fruiting_node(k, l, delay_frt_by_c_stress, stem_n_ratio);
            }
            legacy = crate::LegacyGlobalState::from_globals();
            for m in 0..legacy.num_nodes[[k, l]] as usize {
                fruiting_site(k, l, m, &mut nwfl);
            }
            legacy = crate::LegacyGlobalState::from_globals();
        }
        legacy = crate::LegacyGlobalState::from_globals();
    }
    scratch.nwfl = nwfl;
    with_phenology_scratch_mut(|state| *state = scratch);

    fruiting_sites_abscission();
    leaf_abscission();
}

fn pre_fruiting_node(stem_n_ratio: f64) {
    const MAX_AGE_PRE_FR_NODE: f64 = 66.0;

    let mut legacy = crate::LegacyGlobalState::from_globals();
    let last_index = (legacy.num_pre_fru_nodes - 1) as usize;
    if legacy.age_of_pre_fru_node[last_index] > MAX_AGE_PRE_FR_NODE {
        return;
    }

    for j in 0..legacy.num_pre_fru_nodes as usize {
        legacy.age_of_pre_fru_node[j] += legacy.day_inc;
    }

    if legacy.num_pre_fru_nodes >= 9 {
        legacy.write_to_globals();
        return;
    }

    let mut time_to_next_pre_fru_node = legacy.var_par[31];
    if legacy.num_pre_fru_nodes <= 2 {
        time_to_next_pre_fru_node *= legacy.var_par[32];
    } else if legacy.num_pre_fru_nodes == 3 {
        time_to_next_pre_fru_node *= legacy.var_par[33];
    }

    let last_index = (legacy.num_pre_fru_nodes - 1) as usize;
    if legacy.age_of_pre_fru_node[last_index] >= time_to_next_pre_fru_node {
        legacy.num_pre_fru_nodes += 1;
        let new_index = (legacy.num_pre_fru_nodes - 1) as usize;
        legacy.node_layer_pre_fru[new_index] = ((legacy.plant_height / 5.0) as i32).clamp(0, 19);
        legacy.leaf_area_pre_fru[new_index] = legacy.var_par[34];
        legacy.leaf_weight_pre_fru[new_index] =
            legacy.leaf_area_pre_fru[new_index] * legacy.leaf_weight_area_ratio;
        legacy.total_stem_weight -= legacy.leaf_weight_pre_fru[new_index];
        legacy.leaf_nitrogen += legacy.leaf_weight_pre_fru[new_index] * stem_n_ratio;
        legacy.stem_nitrogen -= legacy.leaf_weight_pre_fru[new_index] * stem_n_ratio;
    }
    legacy.write_to_globals();
}

fn days_to_first_square() -> f64 {
    let legacy = crate::LegacyGlobalState::from_globals();
    with_phenology_scratch_mut(|scratch| {
        let p1 = 34.0;
        let p2 = 132.2;
        let p3 = -7.0;
        let p4 = 0.125;
        let p5 = 0.08;
        let p6 = 0.30;

        if legacy.daynum <= legacy.day_emerge {
            scratch.days_to_1st_avtemp = legacy.avrg_daily_temp;
            scratch.days_to_1st_sumstrs = 0.0;
        }

        scratch.days_to_1st_avtemp = ((legacy.kday - 1) as f64 * scratch.days_to_1st_avtemp
            + legacy.avrg_daily_temp)
            / legacy.kday as f64;
        if scratch.days_to_1st_avtemp > p1 {
            scratch.days_to_1st_avtemp = p1;
        }

        let mut tsq1 = p2 + scratch.days_to_1st_avtemp * (p3 + scratch.days_to_1st_avtemp * p4);
        tsq1 *= legacy.var_par[30];
        scratch.days_to_1st_sumstrs +=
            p5 * (1.0 - legacy.water_stress) + p6 * (1.0 - legacy.n_stress_veg);
        tsq1 -= scratch.days_to_1st_sumstrs;
        tsq1
    })
}

fn create_first_square(stem_n_ratio: f64) {
    let mut legacy = crate::LegacyGlobalState::from_globals();
    legacy.fruiting_code[[0, 0, 0]] = 1;
    legacy.fruit_fraction[[0, 0, 0]] = 1.0;

    legacy.leaf_area_nodes[[0, 0, 0]] = legacy.var_par[34];
    legacy.leaf_weight_nodes[[0, 0, 0]] = legacy.var_par[34] * legacy.leaf_weight_area_ratio;
    legacy.total_stem_weight -= legacy.leaf_weight_nodes[[0, 0, 0]];
    legacy.leaf_nitrogen += legacy.leaf_weight_nodes[[0, 0, 0]] * stem_n_ratio;
    legacy.stem_nitrogen -= legacy.leaf_weight_nodes[[0, 0, 0]] * stem_n_ratio;

    legacy.num_fruit_branches[0] = 1;
    legacy.num_nodes[[0, 0]] = 1;
    legacy.node_layer[[0, 0]] = ((legacy.plant_height / 5.0) as i32).clamp(0, 19);
    legacy.fruit_growth_ratio = 1.0;
    legacy.avrg_node_temper[[0, 0, 0]] = legacy.avrg_daily_temp;

    let cotylwt = 0.20;
    legacy.abscised_leaf_weight += cotylwt;
    legacy.cum_plant_n_loss += cotylwt * legacy.leaf_nitrogen / TotalLeafWeight();
    legacy.leaf_nitrogen -= cotylwt * legacy.leaf_nitrogen / TotalLeafWeight();
    legacy.pix_in_plants -= cotylwt * legacy.pixcon;
    legacy.write_to_globals();
}

fn add_vegetative_branch(delay_veg_by_c_stress: f64, stem_n_ratio: f64, days_to_1st_square: f64) {
    let phen_delay_by_n_stress = read_phenology_scratch().phen_delay_by_n_stress;
    const VPVEGB: [f64; 3] = [13.39, -0.696, 0.012];

    let mut legacy = crate::LegacyGlobalState::from_globals();
    let prev_branch = (legacy.num_veg_branches - 1) as usize;
    let time_to_next_veg_branch = VPVEGB[0]
        + legacy.avrg_node_temper[[prev_branch, 0, 0]]
            * (VPVEGB[1] + legacy.avrg_node_temper[[prev_branch, 0, 0]] * VPVEGB[2]);

    let threshold = time_to_next_veg_branch
        + delay_veg_by_c_stress
        + phen_delay_by_n_stress
        + days_to_1st_square;
    if legacy.age_of_site[[prev_branch, 0, 0]] < threshold {
        return;
    }

    legacy.num_veg_branches += 1;
    if legacy.num_veg_branches > 3 {
        legacy.num_veg_branches = 3;
        legacy.write_to_globals();
        return;
    }

    let new_branch = (legacy.num_veg_branches - 1) as usize;
    legacy.fruit_fraction[[new_branch, 0, 0]] = 1.0;
    legacy.fruiting_code[[new_branch, 0, 0]] = 1;

    legacy.leaf_area_nodes[[new_branch, 0, 0]] = legacy.var_par[34];
    legacy.leaf_weight_nodes[[new_branch, 0, 0]] =
        legacy.var_par[34] * legacy.leaf_weight_area_ratio;
    legacy.leaf_area_main_stem[[new_branch, 0]] = legacy.var_par[34];
    legacy.leaf_weight_main_stem[[new_branch, 0]] =
        legacy.leaf_area_main_stem[[new_branch, 0]] * legacy.leaf_weight_area_ratio;

    legacy.total_stem_weight -= legacy.leaf_weight_nodes[[new_branch, 0, 0]]
        + legacy.leaf_weight_main_stem[[new_branch, 0]];
    let addlfn = (legacy.leaf_weight_nodes[[new_branch, 0, 0]]
        + legacy.leaf_weight_main_stem[[new_branch, 0]])
        * stem_n_ratio;
    legacy.leaf_nitrogen += addlfn;
    legacy.stem_nitrogen -= addlfn;

    legacy.avrg_node_temper[[new_branch, 0, 0]] = legacy.avrg_daily_temp;
    legacy.num_fruit_branches[new_branch] = 1;
    legacy.num_nodes[[new_branch, 0]] = 1;
    legacy.node_layer[[new_branch, 0]] = ((legacy.plant_height / 5.0) as i32).clamp(0, 19);
    legacy.write_to_globals();
}

fn add_fruiting_branch(k: usize, delay_veg_by_c_stress: f64, stem_n_ratio: f64) {
    let phen_delay_by_n_stress = read_phenology_scratch().phen_delay_by_n_stress;
    const VFRTBR: [f64; 8] = [0.8, 0.95, 33.0, 4.461, -0.1912, 0.00265, 1.8, -1.32];

    let mut legacy = crate::LegacyGlobalState::from_globals();
    legacy.delay_new_fru_branch[k] += delay_veg_by_c_stress + VFRTBR[0] * phen_delay_by_n_stress;
    legacy.delay_new_fru_branch[k] += VFRTBR[1] * (1.0 - legacy.water_stress);

    let nbrch = (legacy.num_fruit_branches[k] - 1) as usize;
    let mut tav = legacy.avrg_node_temper[[k, nbrch, 0]];
    if tav > VFRTBR[2] {
        tav = VFRTBR[2];
    }

    let mut time_to_next_fru_branch =
        legacy.var_par[35] + tav * (VFRTBR[3] + tav * (VFRTBR[4] + tav * VFRTBR[5]));
    if k > 0 {
        time_to_next_fru_branch *= VFRTBR[6];
    }
    time_to_next_fru_branch = time_to_next_fru_branch
        * (1.0 + VFRTBR[7] * (1.0 - legacy.density_factor))
        + legacy.delay_new_fru_branch[k];

    if legacy.kday > legacy.kday_adjust
        && legacy.kday <= (legacy.kday_adjust + legacy.num_adjust_days)
        && legacy.nadj[0]
    {
        if legacy.adj_add_ms_nodes_rate == 0.0 {
            time_to_next_fru_branch = 100.0;
        } else {
            time_to_next_fru_branch /= legacy.adj_add_ms_nodes_rate;
        }
    }

    if legacy.age_of_site[[k, nbrch, 0]] < time_to_next_fru_branch {
        return;
    }

    legacy.num_fruit_branches[k] += 1;
    if legacy.num_fruit_branches[k] > 30 {
        legacy.num_fruit_branches[k] = 30;
        legacy.write_to_globals();
        return;
    }

    let newbr = (legacy.num_fruit_branches[k] - 1) as usize;
    legacy.num_nodes[[k, newbr]] = 1;
    legacy.fruit_fraction[[k, newbr, 0]] = 1.0;
    legacy.fruiting_code[[k, newbr, 0]] = 1;

    legacy.leaf_area_nodes[[k, newbr, 0]] = legacy.var_par[34];
    legacy.leaf_weight_nodes[[k, newbr, 0]] = legacy.var_par[34] * legacy.leaf_weight_area_ratio;
    legacy.leaf_area_main_stem[[k, newbr]] = legacy.var_par[34];
    legacy.leaf_weight_main_stem[[k, newbr]] =
        legacy.leaf_area_main_stem[[k, newbr]] * legacy.leaf_weight_area_ratio;

    legacy.total_stem_weight -=
        legacy.leaf_weight_main_stem[[k, newbr]] + legacy.leaf_weight_nodes[[k, newbr, 0]];
    let addlfn = (legacy.leaf_weight_main_stem[[k, newbr]]
        + legacy.leaf_weight_nodes[[k, newbr, 0]])
        * stem_n_ratio;
    legacy.leaf_nitrogen += addlfn;
    legacy.node_layer[[k, newbr]] = ((legacy.plant_height / 5.0) as i32).clamp(0, 19);
    legacy.stem_nitrogen -= addlfn;

    legacy.avrg_node_temper[[k, newbr, 0]] = legacy.avrg_daily_temp;
    legacy.delay_new_fru_branch[k] = 0.0;
    legacy.write_to_globals();
}

fn add_fruiting_node(k: usize, l: usize, delay_frt_by_c_stress: f64, stem_n_ratio: f64) {
    let phen_delay_by_n_stress = read_phenology_scratch().phen_delay_by_n_stress;
    const VFRTNOD: [f64; 6] = [1.32, 0.90, 33.0, 7.6725, -0.3297, 0.004657];

    let mut legacy = crate::LegacyGlobalState::from_globals();
    legacy.delay_new_node[[k, l]] +=
        (delay_frt_by_c_stress + VFRTNOD[0] * phen_delay_by_n_stress) / legacy.pixdn;
    legacy.delay_new_node[[k, l]] += VFRTNOD[1] * (1.0 - legacy.water_stress);

    let nnid = (legacy.num_nodes[[k, l]] - 1) as usize;
    let mut tav = legacy.avrg_node_temper[[k, l, nnid]];
    if tav > VFRTNOD[2] {
        tav = VFRTNOD[2];
    }

    let mut time_to_next_fru_node =
        legacy.var_par[36] + tav * (VFRTNOD[3] + tav * (VFRTNOD[4] + tav * VFRTNOD[5]));
    time_to_next_fru_node = time_to_next_fru_node
        * (1.0 + legacy.var_par[37] * (1.0 - legacy.density_factor))
        + legacy.delay_new_node[[k, l]];

    if legacy.kday > legacy.kday_adjust
        && legacy.kday <= (legacy.kday_adjust + legacy.num_adjust_days)
        && legacy.nadj[2]
    {
        if legacy.adj_add_sites_rate == 0.0 {
            time_to_next_fru_node = 100.0;
        } else {
            time_to_next_fru_node /= legacy.adj_add_sites_rate;
        }
    }

    if legacy.age_of_site[[k, l, nnid]] < time_to_next_fru_node {
        return;
    }

    legacy.num_nodes[[k, l]] += 1;
    if legacy.num_nodes[[k, l]] > 5 {
        legacy.num_nodes[[k, l]] = 5;
        legacy.write_to_globals();
        return;
    }

    let newnod = nnid + 1;
    legacy.fruit_fraction[[k, l, newnod]] = 1.0;
    legacy.fruiting_code[[k, l, newnod]] = 1;
    legacy.leaf_area_nodes[[k, l, newnod]] = legacy.var_par[34];
    legacy.leaf_weight_nodes[[k, l, newnod]] = legacy.var_par[34] * legacy.leaf_weight_area_ratio;
    legacy.total_stem_weight -= legacy.leaf_weight_nodes[[k, l, newnod]];
    legacy.leaf_nitrogen += legacy.leaf_weight_nodes[[k, l, newnod]] * stem_n_ratio;
    legacy.stem_nitrogen -= legacy.leaf_weight_nodes[[k, l, newnod]] * stem_n_ratio;

    legacy.avrg_node_temper[[k, l, newnod]] = legacy.avrg_daily_temp;
    legacy.delay_new_node[[k, l]] = 0.0;
    legacy.write_to_globals();
}

fn fruiting_site(k: usize, l: usize, m: usize, node_recent_white_flower: &mut i32) {
    let mut scratch = read_phenology_scratch();
    const VFRSITE: [f64; 15] = [
        0.60, 0.40, 12.25, 0.40, 33.0, 0.20, 0.04, 0.45, 26.10, 9.0, 0.10, 3.0, 1.129, 0.043, 0.26,
    ];

    let mut legacy = crate::LegacyGlobalState::from_globals();
    if legacy.fruiting_code[[k, l, m]] <= 0 {
        scratch.boltmp[k][l][m] = 0.0;
        with_phenology_scratch_mut(|state| *state = scratch);
        return;
    }

    if legacy.num_fruit_sites <= 0 {
        scratch.fib_length = 0.0;
        scratch.fib_strength = 0.0;
    }
    legacy.num_fruit_sites += 1;

    let agefac =
        (1.0 - legacy.water_stress) * VFRSITE[0] + (1.0 - legacy.n_stress_veg) * VFRSITE[1];
    legacy.leaf_age[[k, l, m]] += legacy.day_inc + agefac;
    if legacy.day_first_def > 0 && legacy.daynum > legacy.day_first_def {
        legacy.leaf_age[[k, l, m]] += legacy.var_par[38];
    }

    if legacy.fruiting_code[[k, l, m]] >= 3 && legacy.fruiting_code[[k, l, m]] <= 6 {
        legacy.write_to_globals();
        with_phenology_scratch_mut(|state| *state = scratch);
        return;
    }

    let tmin = GetFromClim(CLIMATE_METRIC_TMIN, legacy.daynum);
    let tmax = GetFromClim(CLIMATE_METRIC_TMAX, legacy.daynum);
    let mut ageinc = legacy.day_inc;
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

    legacy.avrg_node_temper[[k, l, m]] = (legacy.avrg_node_temper[[k, l, m]]
        * legacy.age_of_site[[k, l, m]]
        + legacy.avrg_daily_temp * ageinc)
        / (legacy.age_of_site[[k, l, m]] + ageinc);
    legacy.age_of_site[[k, l, m]] += ageinc;

    if legacy.fruiting_code[[k, l, m]] == 1 {
        if legacy.age_of_site[[k, l, m]] >= VFRSITE[8] {
            scratch.boltmp[k][l][m] = legacy.avrg_daily_temp;
            legacy.age_of_boll[[k, l, m]] = legacy.day_inc;
            legacy.fruiting_code[[k, l, m]] = 7;
            legacy.write_to_globals();
            new_boll_formation(k, l, m);
            legacy = crate::LegacyGlobalState::from_globals();
            if legacy.cotton_weight_green_bolls > 0.0 && legacy.first_bloom <= 1 {
                legacy.first_bloom = legacy.daynum;
            }
            if k == 0 && m == 0 && l as i32 > *node_recent_white_flower {
                *node_recent_white_flower = l as i32;
            }
        }
        legacy.write_to_globals();
        with_phenology_scratch_mut(|state| *state = scratch);
        return;
    }

    if legacy.boll_weight[[k, l, m]] > 0.0 {
        let dum = if legacy.leaf_area_index <= VFRSITE[11] && legacy.kday > 100 {
            VFRSITE[12] - VFRSITE[13] * legacy.leaf_area_index
        } else {
            1.0
        };
        let dagebol = legacy.day_inc * dum
            + VFRSITE[14] * (1.0 - legacy.water_stress)
            + VFRSITE[10] * (1.0 - legacy.n_stress_fruiting);
        scratch.boltmp[k][l][m] = (scratch.boltmp[k][l][m] * legacy.age_of_boll[[k, l, m]]
            + legacy.avrg_daily_temp * dagebol)
            / (legacy.age_of_boll[[k, l, m]] + dagebol);
        legacy.age_of_boll[[k, l, m]] += dagebol;
    }

    if legacy.fruiting_code[[k, l, m]] == 7 {
        if legacy.age_of_boll[[k, l, m]] >= VFRSITE[9] {
            legacy.fruiting_code[[k, l, m]] = 2;
        }
        legacy.write_to_globals();
        with_phenology_scratch_mut(|state| *state = scratch);
        return;
    }

    if legacy.fruiting_code[[k, l, m]] == 2 {
        legacy.write_to_globals();
        boll_opening(k, l, m, scratch.boltmp[k][l][m]);
    } else {
        legacy.write_to_globals();
    }
    with_phenology_scratch_mut(|state| *state = scratch);
}

fn new_boll_formation(k: usize, l: usize, m: usize) {
    const SEED_RATIO: f64 = 0.64;
    const VNEWBOLL: [f64; 2] = [0.31, 0.02];

    let mut legacy = crate::LegacyGlobalState::from_globals();
    if !legacy.b_pollin_switch {
        legacy.fruiting_code[[k, l, m]] = 6;
        legacy.fruit_fraction[[k, l, m]] = 0.0;
        legacy.bloom_weight_loss += legacy.square_weight[[k, l, m]];
        legacy.square_weight[[k, l, m]] = 0.0;
        legacy.write_to_globals();
        return;
    }

    let bolinit = VNEWBOLL[0] * legacy.square_weight[[k, l, m]];
    legacy.boll_weight[[k, l, m]] = 0.2 * bolinit;
    legacy.burr_weight[[k, l, m]] = bolinit - legacy.boll_weight[[k, l, m]];
    legacy.bloom_weight_loss += legacy.square_weight[[k, l, m]] - bolinit;

    let mut sqr1n = legacy.square_n_conc * legacy.square_weight[[k, l, m]];
    legacy.square_nitrogen -= sqr1n;
    legacy.cum_plant_n_loss += sqr1n * (1.0 - VNEWBOLL[0]);
    sqr1n *= VNEWBOLL[0];

    let mut seed1n = legacy.boll_weight[[k, l, m]] * SEED_RATIO * VNEWBOLL[1];
    if seed1n > sqr1n {
        seed1n = sqr1n;
    }
    legacy.seed_nitrogen += seed1n;
    legacy.burr_nitrogen += sqr1n - seed1n;

    legacy.cotton_weight_green_bolls += legacy.boll_weight[[k, l, m]];
    legacy.burr_weight_green_bolls += legacy.burr_weight[[k, l, m]];
    legacy.total_square_weight -= legacy.square_weight[[k, l, m]];
    legacy.square_weight[[k, l, m]] = 0.0;
    legacy.write_to_globals();
}

fn boll_opening(k: usize, l: usize, m: usize, tmpboll: f64) {
    let ddpar1 = 1.0;
    let ddpar2 = 0.8;
    let vboldhs: [f64; 11] = [
        30.0, 41.189, -1.6057, 0.020743, 70.0, 0.994, 56.603, -2.921, 0.059, 1.219, 0.0065,
    ];

    let mut legacy = crate::LegacyGlobalState::from_globals();
    let mut atn = tmpboll;
    if atn > vboldhs[0] {
        atn = vboldhs[0];
    }

    let mut dehiss =
        legacy.var_par[39] + atn * (vboldhs[1] + atn * (vboldhs[2] + atn * vboldhs[3]));
    dehiss *= legacy.var_par[40];
    if dehiss > vboldhs[4] {
        dehiss = vboldhs[4];
    }

    if legacy.day_first_def > 0 && legacy.daynum > legacy.day_first_def {
        dehiss *= vboldhs[5].powf((legacy.daynum - legacy.day_first_def) as f64);
    }

    if legacy.leaf_area_index < ddpar1 {
        let mut fdhslai = ddpar2 + legacy.leaf_area_index * (1.0 - ddpar2) / ddpar1;
        fdhslai = fdhslai.clamp(0.0, 1.0);
        dehiss *= fdhslai;
    }

    if legacy.age_of_boll[[k, l, m]] < dehiss {
        return;
    }

    legacy.fruiting_code[[k, l, m]] = 3;
    legacy.cotton_weight_open_bolls += legacy.boll_weight[[k, l, m]];
    legacy.burr_weight_open_bolls += legacy.burr_weight[[k, l, m]];
    legacy.cotton_weight_green_bolls -= legacy.boll_weight[[k, l, m]];
    legacy.burr_weight_green_bolls -= legacy.burr_weight[[k, l, m]];

    legacy.ginp = (legacy.var_par[41] - legacy.var_par[42] * atn) / 100.0;
    legacy.gintot = (legacy.gintot * legacy.num_open_bolls
        + legacy.ginp * legacy.fruit_fraction[[k, l, m]])
        / (legacy.num_open_bolls + legacy.fruit_fraction[[k, l, m]]);
    legacy.lint_yield +=
        legacy.ginp * legacy.boll_weight[[k, l, m]] * legacy.plant_population * 0.001;

    let fsx = vboldhs[6] + atn * (vboldhs[7] + vboldhs[8] * atn);
    let flx = vboldhs[9] - vboldhs[10] * atn;
    with_phenology_scratch_mut(|scratch| {
        scratch.fib_strength = (scratch.fib_strength * legacy.num_open_bolls
            + fsx * legacy.fruit_fraction[[k, l, m]])
            / (legacy.num_open_bolls + legacy.fruit_fraction[[k, l, m]]);
        scratch.fib_length = (scratch.fib_length * legacy.num_open_bolls
            + flx * legacy.fruit_fraction[[k, l, m]])
            / (legacy.num_open_bolls + legacy.fruit_fraction[[k, l, m]]);
    });

    legacy.num_open_bolls += legacy.fruit_fraction[[k, l, m]];
    legacy.write_to_globals();
}
