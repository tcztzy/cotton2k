use crate::model_state::for_each_row;
use crate::LegacyGlobalState;
use ndarray::{s, Array1, Array2, ArrayView1, ArrayViewMut1};
use std::sync::{LazyLock, RwLock};

static SOIL_N_DEPTH: LazyLock<RwLock<[f64; 40]>> = LazyLock::new(|| RwLock::new([0.0; 40]));

struct SoilNitrogenState<'a> {
    daynum: i32,
    day_start: i32,
    nk: usize,
    nl: usize,
    soil_n_depth: &'a mut [f64; 40],
    dl: &'a Array1<f64>,
    wk: &'a Array1<f64>,
    thetar: &'a Array1<f64>,
    thts: &'a Array1<f64>,
    bulk_density: &'a Array1<f64>,
    field_capacity: &'a Array1<f64>,
    soil_horizon_num: &'a Array1<i32>,
    fresh_organic_matter: &'a mut Array2<f64>,
    fresh_organic_nitrogen: &'a mut Array2<f64>,
    humus_nitrogen: &'a mut Array2<f64>,
    humus_organic_matter: &'a mut Array2<f64>,
    soil_temp_daily_avrg: &'a Array2<f64>,
    vol_nh4_n_content: &'a mut Array2<f64>,
    vol_no3_n_content: &'a mut Array2<f64>,
    vol_urea_n_content: &'a mut Array2<f64>,
    vol_water_content: &'a Array2<f64>,
    mineralized_organic_n: &'a mut f64,
    soil_nitrogen_loss: &'a mut f64,
}

struct SoilNitrogenRowView<'a> {
    fresh_organic_matter: ArrayViewMut1<'a, f64>,
    fresh_organic_nitrogen: ArrayViewMut1<'a, f64>,
    humus_nitrogen: ArrayViewMut1<'a, f64>,
    humus_organic_matter: ArrayViewMut1<'a, f64>,
    soil_temp_daily_avrg: ArrayView1<'a, f64>,
    vol_nh4_n_content: ArrayViewMut1<'a, f64>,
    vol_no3_n_content: ArrayViewMut1<'a, f64>,
    vol_urea_n_content: ArrayViewMut1<'a, f64>,
    vol_water_content: ArrayView1<'a, f64>,
    wk: &'a Array1<f64>,
    dl: f64,
    thetar: f64,
    thts: f64,
    field_capacity: f64,
    bulk_density: f64,
    soil_depth: f64,
}

pub fn soil_nitrogen() {
    let mut soil_n_depth = SOIL_N_DEPTH
        .write()
        .expect("soil nitrogen depth state lock should not be poisoned");
    let mut legacy = LegacyGlobalState::from_globals();
    let mut state = SoilNitrogenState {
        daynum: legacy.daynum,
        day_start: legacy.day_start,
        nk: legacy.nk as usize,
        nl: legacy.nl as usize,
        soil_n_depth: &mut soil_n_depth,
        dl: &legacy.dl,
        wk: &legacy.wk,
        thetar: &legacy.thetar,
        thts: &legacy.thts,
        bulk_density: &legacy.bulk_density,
        field_capacity: &legacy.field_capacity,
        soil_horizon_num: &legacy.soil_horizon_num,
        fresh_organic_matter: &mut legacy.fresh_organic_matter,
        fresh_organic_nitrogen: &mut legacy.fresh_organic_nitrogen,
        humus_nitrogen: &mut legacy.humus_nitrogen,
        humus_organic_matter: &mut legacy.humus_organic_matter,
        soil_temp_daily_avrg: &legacy.soil_temp_daily_avrg,
        vol_nh4_n_content: &mut legacy.vol_nh4_n_content,
        vol_no3_n_content: &mut legacy.vol_no3_n_content,
        vol_urea_n_content: &mut legacy.vol_urea_n_content,
        vol_water_content: &legacy.vol_water_content,
        mineralized_organic_n: &mut legacy.mineralized_organic_n,
        soil_nitrogen_loss: &mut legacy.soil_nitrogen_loss,
    };
    soil_nitrogen_impl(&mut state);
    legacy.write_to_globals();
}

fn soil_nitrogen_impl(state: &mut SoilNitrogenState<'_>) {
    let initialize_organic_n = state.daynum <= state.day_start;
    if initialize_organic_n {
        initialize_soil_n_depth(state.nl, state.dl, state.soil_n_depth);
    }

    let mut mineralized_organic_n_delta = 0.0;
    let mut soil_nitrogen_loss_delta = 0.0;

    for layer in 0..state.nl {
        let horizon = state.soil_horizon_num[layer] as usize;
        let mut row = SoilNitrogenRowView {
            fresh_organic_matter: state.fresh_organic_matter.slice_mut(s![layer, ..]),
            fresh_organic_nitrogen: state.fresh_organic_nitrogen.slice_mut(s![layer, ..]),
            humus_nitrogen: state.humus_nitrogen.slice_mut(s![layer, ..]),
            humus_organic_matter: state.humus_organic_matter.slice_mut(s![layer, ..]),
            soil_temp_daily_avrg: state.soil_temp_daily_avrg.slice(s![layer, ..]),
            vol_nh4_n_content: state.vol_nh4_n_content.slice_mut(s![layer, ..]),
            vol_no3_n_content: state.vol_no3_n_content.slice_mut(s![layer, ..]),
            vol_urea_n_content: state.vol_urea_n_content.slice_mut(s![layer, ..]),
            vol_water_content: state.vol_water_content.slice(s![layer, ..]),
            wk: state.wk,
            dl: state.dl[layer],
            thetar: state.thetar[layer],
            thts: state.thts[layer],
            field_capacity: state.field_capacity[layer],
            bulk_density: state.bulk_density[horizon],
            soil_depth: state.soil_n_depth[layer],
        };

        for column in 0..state.nk {
            let urea = row.vol_urea_n_content[column];
            let water = row.vol_water_content[column];
            let temp = row.soil_temp_daily_avrg[column];
            if urea > 0.0 {
                urea_hydrolysis(&mut row, column);
            }
            mineralized_organic_n_delta +=
                mineralize_nitrogen(&mut row, column, initialize_organic_n);
            if row.vol_nh4_n_content[column] > 0.00001 {
                nitrification(&mut row, column);
            }

            if should_denitrify(
                row.vol_no3_n_content[column],
                water,
                row.field_capacity,
                temp,
            ) {
                soil_nitrogen_loss_delta += denitrification(&mut row, column);
            }
        }
    }

    *state.mineralized_organic_n += mineralized_organic_n_delta;
    *state.soil_nitrogen_loss += soil_nitrogen_loss_delta;
}

fn initialize_soil_n_depth(
    n_layers: usize,
    layer_depths: &Array1<f64>,
    soil_n_depth: &mut [f64; 40],
) {
    let mut sum_depth = 0.0;
    for (layer, &layer_depth) in layer_depths.iter().take(n_layers).enumerate() {
        sum_depth += layer_depth;
        soil_n_depth[layer] = sum_depth;
    }
}

fn should_denitrify(no3: f64, water: f64, field_capacity: f64, soil_temp_k: f64) -> bool {
    no3 > 0.001 && water > field_capacity && soil_temp_k >= 278.161
}

fn urea_hydrolysis(row: &mut SoilNitrogenRowView<'_>, column: usize) {
    const AK0: f64 = 0.25;
    const CAK1: f64 = 0.3416;
    const CAK2: f64 = 0.0776;
    const STF1: f64 = 40.0;
    const STF2: f64 = 0.20;
    const SWF1: f64 = 0.20;

    let oc = 0.4 * (row.fresh_organic_matter[column] + row.humus_organic_matter[column]) * 0.1
        / row.bulk_density;

    let mut ak = CAK1 + CAK2 * oc;
    if ak < AK0 {
        ak = AK0;
    }

    let mut swf = soil_water_effect(row, column, 0.5) + SWF1;
    swf = swf.clamp(0.0, 1.0);

    let mut stf = (row.soil_temp_daily_avrg[column] - 273.161) / STF1 + STF2;
    stf = stf.clamp(0.0, 1.0);

    let mut hydrur = ak * swf * stf * row.vol_urea_n_content[column];
    if hydrur > row.vol_urea_n_content[column] {
        hydrur = row.vol_urea_n_content[column];
    }
    row.vol_urea_n_content[column] -= hydrur;
    row.vol_nh4_n_content[column] += hydrur;
}

fn soil_water_effect(row: &SoilNitrogenRowView<'_>, column: usize, xx: f64) -> f64 {
    let mut wf = if row.vol_water_content[column] <= row.field_capacity {
        (row.vol_water_content[column] - row.thetar) / (row.field_capacity - row.thetar)
    } else {
        1.0 - xx * (row.vol_water_content[column] - row.field_capacity)
            / (row.thts - row.field_capacity)
    };

    if wf < 0.0 {
        wf = 0.0;
    }
    wf
}

fn mineralize_nitrogen(
    row: &mut SoilNitrogenRowView<'_>,
    column: usize,
    initialize_organic_n: bool,
) -> f64 {
    const CN_FRESH: f64 = 25.0;
    const CN_HUM: f64 = 10.0;
    const CN_MAX: f64 = 13.0;
    const CPAR_CN_RF: f64 = 0.693;
    const CPAR_HUMUS_N: f64 = 0.20;
    const CPAR_MIN_NH4: f64 = 0.00025;
    const DECAY_RATE_FRESH: f64 = 0.03;
    const DECAY_RATE_HUMUS: f64 = 0.000083;

    if initialize_organic_n {
        row.fresh_organic_nitrogen[column] = row.fresh_organic_matter[column] * 0.4 / CN_FRESH;
        row.humus_nitrogen[column] = row.humus_organic_matter[column] * 0.4 / CN_HUM;
    }

    if row.fresh_organic_matter[column] <= 0.0 && row.humus_organic_matter[column] <= 0.0 {
        return 0.0;
    }

    let mut cn_ratio_effect = 1.0;
    let total_soil_n = row.fresh_organic_nitrogen[column]
        + row.vol_no3_n_content[column]
        + row.vol_nh4_n_content[column];
    if total_soil_n > 0.0 {
        let cn_ratio = row.fresh_organic_matter[column] * 0.4 / total_soil_n;
        if cn_ratio >= 1000.0 {
            cn_ratio_effect = 0.0;
        } else if cn_ratio > CN_MAX {
            cn_ratio_effect = (-CPAR_CN_RF * (cn_ratio - CN_MAX) / CN_MAX).exp();
        }
    }

    let wf = soil_water_effect(row, column, 0.5);
    let tfac = soil_temperature_effect(row.soil_temp_daily_avrg[column] - 273.161);

    let gross_release_n;
    let immobilization_rate_n;
    if row.fresh_organic_matter[column] > 0.00001 {
        let g1 = tfac * wf * cn_ratio_effect * DECAY_RATE_FRESH;
        let gross_release_dw = g1 * row.fresh_organic_matter[column];
        gross_release_n = g1 * row.fresh_organic_nitrogen[column];

        const CPAR_N_REQ: f64 = 0.0165;
        let mut immobilization = gross_release_dw
            * (CPAR_N_REQ - row.fresh_organic_nitrogen[column] / row.fresh_organic_matter[column]);
        let rnac1 =
            row.vol_nh4_n_content[column] + row.vol_no3_n_content[column] - 2.0 * CPAR_MIN_NH4;
        if immobilization > rnac1 {
            immobilization = rnac1;
        }
        if immobilization < 0.0 {
            immobilization = 0.0;
        }

        row.fresh_organic_matter[column] -= gross_release_dw;
        row.fresh_organic_nitrogen[column] += immobilization - gross_release_n;
        immobilization_rate_n = immobilization;
    } else {
        gross_release_n = 0.0;
        immobilization_rate_n = 0.0;
    }

    let rhmin = row.humus_nitrogen[column] * DECAY_RATE_HUMUS * tfac * wf;
    row.humus_nitrogen[column] -= rhmin + CPAR_HUMUS_N * gross_release_n;
    row.humus_organic_matter[column] -=
        CN_HUM * rhmin / 0.4 + CPAR_HUMUS_N * CN_FRESH * gross_release_n / 0.4;

    let net_n_released = (1.0 - CPAR_HUMUS_N) * gross_release_n + rhmin - immobilization_rate_n;
    let mut mineralized_organic_n_delta = 0.0;
    if net_n_released > 0.0 {
        row.vol_nh4_n_content[column] += net_n_released;
        mineralized_organic_n_delta += net_n_released * row.dl * row.wk[column];
    } else {
        let mut nnom1 = 0.0;
        if row.vol_nh4_n_content[column] > CPAR_MIN_NH4 {
            let addvnc = if net_n_released.abs() < (row.vol_nh4_n_content[column] - CPAR_MIN_NH4) {
                -net_n_released
            } else {
                row.vol_nh4_n_content[column] - CPAR_MIN_NH4
            };
            row.vol_nh4_n_content[column] -= addvnc;
            mineralized_organic_n_delta -= addvnc * row.dl * row.wk[column];
            row.fresh_organic_nitrogen[column] += addvnc;
            nnom1 = net_n_released + addvnc;
        }

        if nnom1 < 0.0 && row.vol_no3_n_content[column] > CPAR_MIN_NH4 {
            let addvnc = if nnom1.abs() < (row.vol_no3_n_content[column] - CPAR_MIN_NH4) {
                -nnom1
            } else {
                row.vol_no3_n_content[column] - CPAR_MIN_NH4
            };
            row.vol_no3_n_content[column] -= addvnc;
            row.fresh_organic_nitrogen[column] += addvnc;
            mineralized_organic_n_delta -= addvnc * row.dl * row.wk[column];
        }
    }

    mineralized_organic_n_delta
}

fn soil_temperature_effect(tt: f64) -> f64 {
    let mut tfm = 0.010645 * (0.12979 * tt).exp();
    tfm = tfm.clamp(0.0, 2.0);
    tfm
}

fn nitrification(row: &mut SoilNitrogenRowView<'_>, column: usize) {
    const CPAR_DEPTH: f64 = 0.45;
    const CPAR_NIT1: f64 = 24.635;
    const CPAR_NIT2: f64 = 8227.0;
    const CPAR_SANC: f64 = 204.0;

    let sanc = if row.vol_nh4_n_content[column] < 0.1 {
        1.0 - (-CPAR_SANC * row.vol_nh4_n_content[column]).exp()
    } else {
        1.0
    };

    let con1 = (CPAR_NIT1 - CPAR_NIT2 / row.soil_temp_daily_avrg[column]).exp();
    let mut ratenit = 1.0 - (-con1).exp();
    let mut tff = (row.soil_depth - 30.0) / 30.0;
    if tff < 0.0 {
        tff = 0.0;
    }

    ratenit = ratenit * sanc * soil_water_effect(row, column, 1.0) * CPAR_DEPTH.powf(tff);
    ratenit = ratenit.clamp(0.0, 0.10);

    let dnit = ratenit * row.vol_nh4_n_content[column];
    row.vol_nh4_n_content[column] -= dnit;
    row.vol_no3_n_content[column] += dnit;
}

fn denitrification(row: &mut SoilNitrogenRowView<'_>, column: usize) -> f64 {
    const CPAR_01: f64 = 24.5;
    const CPAR_02: f64 = 3.1;
    const CPAR_DENIT: f64 = 0.00006;
    const CPAR_FT: f64 = 0.046;
    const CPAR_HUM: f64 = 0.58;
    const VNO3_MIN: f64 = 0.00025;

    let soilc = CPAR_HUM * row.humus_organic_matter[column];
    let cw = CPAR_01 + CPAR_02 * soilc;

    let mut fw =
        (row.vol_water_content[column] - row.field_capacity) / (row.thts - row.field_capacity);
    if fw < 0.0 {
        fw = 0.0;
    }

    let mut ft = 0.1 * (CPAR_FT * (row.soil_temp_daily_avrg[column] - 273.161)).exp();
    if ft > 1.0 {
        ft = 1.0;
    }

    let mut dnrate = CPAR_DENIT * cw * row.vol_no3_n_content[column] * fw * ft;
    if dnrate > row.vol_no3_n_content[column] - VNO3_MIN {
        dnrate = row.vol_no3_n_content[column] - VNO3_MIN;
    }
    if dnrate < 0.0 {
        dnrate = 0.0;
    }

    row.vol_no3_n_content[column] -= dnrate;
    dnrate * row.dl * row.wk[column]
}

pub fn soil_nitrogen_bal() {
    let legacy = LegacyGlobalState::from_globals();
    let balsn =
        legacy.soil_nitrogen_at_start + legacy.cum_fertilizer_n + legacy.mineralized_organic_n
            - legacy.cum_nitrogen_uptake
            - legacy.total_soil_nitrogen
            - legacy.soil_nitrogen_loss;
    let _ = balsn;
}

pub fn soil_nitrogen_average() {
    let legacy = LegacyGlobalState::from_globals();
    let nk = legacy.nk as usize;
    let mut avno = [0.0; 4];
    let mut avnh = [0.0; 4];

    let band_index = |layer: usize| -> Option<usize> {
        match layer {
            0..=7 => Some(0),
            8..=13 => Some(1),
            14..=19 => Some(2),
            20..=25 => Some(3),
            _ => None,
        }
    };

    for_each_row(&legacy.vol_no3_n_content, 26, |layer, row| {
        if let Some(band) = band_index(layer) {
            let row_sum: f64 = row.iter().take(nk).sum();
            avno[band] += row_sum * legacy.dl[layer];
        }
    });
    for_each_row(&legacy.vol_nh4_n_content, 26, |layer, row| {
        if let Some(band) = band_index(layer) {
            let row_sum: f64 = row.iter().take(nk).sum();
            avnh[band] += row_sum * legacy.dl[layer];
        }
    });

    let mut avno30 = avno[0];
    let mut avno60 = avno[1];
    let mut avno90 = avno[2];
    let mut avno120 = avno[3];
    let mut avnh30 = avnh[0];
    let mut avnh60 = avnh[1];
    let mut avnh90 = avnh[2];
    let mut avnh120 = avnh[3];

    avno30 = 1000.0 * avno30 / (30.0 * legacy.nk as f64);
    avnh30 = 1000.0 * avnh30 / (30.0 * legacy.nk as f64);
    avno60 = 1000.0 * avno60 / (30.0 * legacy.nk as f64);
    avnh60 = 1000.0 * avnh60 / (30.0 * legacy.nk as f64);
    avno90 = 1000.0 * avno90 / (30.0 * legacy.nk as f64);
    avnh90 = 1000.0 * avnh90 / (30.0 * legacy.nk as f64);
    avno120 = 1000.0 * avno120 / (30.0 * legacy.nk as f64);
    avnh120 = 1000.0 * avnh120 / (30.0 * legacy.nk as f64);

    let _ = (
        avno30, avnh30, avno60, avnh60, avno90, avnh90, avno120, avnh120,
    );
}
