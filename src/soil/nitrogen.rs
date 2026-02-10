use crate::{
    dl, nk, nl, thetar, thts, wk, BulkDensity, CumFertilizerN, CumNitrogenUptake, DayStart, Daynum,
    FieldCapacity, FreshOrganicMatter, FreshOrganicNitrogen, HumusNitrogen, HumusOrganicMatter,
    MineralizedOrganicN, SoilHorizonNum, SoilNitrogenAtStart, SoilNitrogenLoss, SoilTempDailyAvrg,
    TotalSoilNitrogen, VolNh4NContent, VolNo3NContent, VolUreaNContent, VolWaterContent,
};

static mut SOIL_N_DEPTH: [f64; 40] = [0.0; 40];

struct SoilNitrogenState<'a> {
    daynum: i32,
    day_start: i32,
    nk: usize,
    nl: usize,
    soil_n_depth: &'a mut [f64; 40],
    dl: &'a [f64; 40],
    wk: &'a [f64; 20],
    thetar: &'a [f64; 40],
    thts: &'a [f64; 40],
    bulk_density: &'a [f64; 9],
    field_capacity: &'a [f64; 40],
    soil_horizon_num: &'a [i32; 40],
    fresh_organic_matter: &'a mut [[f64; 20]; 40],
    fresh_organic_nitrogen: &'a mut [[f64; 20]; 40],
    humus_nitrogen: &'a mut [[f64; 20]; 40],
    humus_organic_matter: &'a mut [[f64; 20]; 40],
    soil_temp_daily_avrg: &'a [[f64; 20]; 40],
    vol_nh4_n_content: &'a mut [[f64; 20]; 40],
    vol_no3_n_content: &'a mut [[f64; 20]; 40],
    vol_urea_n_content: &'a mut [[f64; 20]; 40],
    vol_water_content: &'a [[f64; 20]; 40],
    mineralized_organic_n: &'a mut f64,
    soil_nitrogen_loss: &'a mut f64,
}

#[allow(static_mut_refs)]
pub fn soil_nitrogen() {
    unsafe {
        let mut state = SoilNitrogenState {
            daynum: Daynum,
            day_start: DayStart,
            nk: nk as usize,
            nl: nl as usize,
            soil_n_depth: &mut SOIL_N_DEPTH,
            dl: &dl,
            wk: &wk,
            thetar: &thetar,
            thts: &thts,
            bulk_density: &BulkDensity,
            field_capacity: &FieldCapacity,
            soil_horizon_num: &SoilHorizonNum,
            fresh_organic_matter: &mut FreshOrganicMatter,
            fresh_organic_nitrogen: &mut FreshOrganicNitrogen,
            humus_nitrogen: &mut HumusNitrogen,
            humus_organic_matter: &mut HumusOrganicMatter,
            soil_temp_daily_avrg: &SoilTempDailyAvrg,
            vol_nh4_n_content: &mut VolNh4NContent,
            vol_no3_n_content: &mut VolNo3NContent,
            vol_urea_n_content: &mut VolUreaNContent,
            vol_water_content: &VolWaterContent,
            mineralized_organic_n: &mut MineralizedOrganicN,
            soil_nitrogen_loss: &mut SoilNitrogenLoss,
        };
        soil_nitrogen_impl(&mut state);
    }
}

fn soil_nitrogen_impl(state: &mut SoilNitrogenState<'_>) {
    if state.daynum <= state.day_start {
        let mut sumdl = 0.0;
        for layer in 0..state.nl {
            sumdl += state.dl[layer];
            state.soil_n_depth[layer] = sumdl;
        }
    }

    for layer in 0..state.nl {
        for column in 0..state.nk {
            if state.vol_urea_n_content[layer][column] > 0.0 {
                urea_hydrolysis(state, layer, column);
            }
            mineralize_nitrogen(state, layer, column);
            if state.vol_nh4_n_content[layer][column] > 0.00001 {
                nitrification(state, layer, column, state.soil_n_depth[layer]);
            }

            if state.vol_no3_n_content[layer][column] > 0.001
                && state.vol_water_content[layer][column] > state.field_capacity[layer]
                && state.soil_temp_daily_avrg[layer][column] >= 278.161
            {
                denitrification(state, layer, column);
            }
        }
    }
}

fn urea_hydrolysis(state: &mut SoilNitrogenState<'_>, layer: usize, column: usize) {
    const AK0: f64 = 0.25;
    const CAK1: f64 = 0.3416;
    const CAK2: f64 = 0.0776;
    const STF1: f64 = 40.0;
    const STF2: f64 = 0.20;
    const SWF1: f64 = 0.20;

    let horizon = state.soil_horizon_num[layer] as usize;
    let oc = 0.4
        * (state.fresh_organic_matter[layer][column] + state.humus_organic_matter[layer][column])
        * 0.1
        / state.bulk_density[horizon];

    let mut ak = CAK1 + CAK2 * oc;
    if ak < AK0 {
        ak = AK0;
    }

    let mut swf = soil_water_effect(state, layer, column, 0.5) + SWF1;
    swf = swf.clamp(0.0, 1.0);

    let mut stf = (state.soil_temp_daily_avrg[layer][column] - 273.161) / STF1 + STF2;
    stf = stf.clamp(0.0, 1.0);

    let mut hydrur = ak * swf * stf * state.vol_urea_n_content[layer][column];
    if hydrur > state.vol_urea_n_content[layer][column] {
        hydrur = state.vol_urea_n_content[layer][column];
    }
    state.vol_urea_n_content[layer][column] -= hydrur;
    state.vol_nh4_n_content[layer][column] += hydrur;
}

fn soil_water_effect(state: &SoilNitrogenState<'_>, layer: usize, column: usize, xx: f64) -> f64 {
    let mut wf = if state.vol_water_content[layer][column] <= state.field_capacity[layer] {
        (state.vol_water_content[layer][column] - state.thetar[layer])
            / (state.field_capacity[layer] - state.thetar[layer])
    } else {
        1.0 - xx * (state.vol_water_content[layer][column] - state.field_capacity[layer])
            / (state.thts[layer] - state.field_capacity[layer])
    };

    if wf < 0.0 {
        wf = 0.0;
    }
    wf
}

fn mineralize_nitrogen(state: &mut SoilNitrogenState<'_>, layer: usize, column: usize) {
    const CN_FRESH: f64 = 25.0;
    const CN_HUM: f64 = 10.0;
    const CN_MAX: f64 = 13.0;
    const CPAR_CN_RF: f64 = 0.693;
    const CPAR_HUMUS_N: f64 = 0.20;
    const CPAR_MIN_NH4: f64 = 0.00025;
    const DECAY_RATE_FRESH: f64 = 0.03;
    const DECAY_RATE_HUMUS: f64 = 0.000083;

    if state.daynum <= state.day_start {
        state.fresh_organic_nitrogen[layer][column] =
            state.fresh_organic_matter[layer][column] * 0.4 / CN_FRESH;
        state.humus_nitrogen[layer][column] =
            state.humus_organic_matter[layer][column] * 0.4 / CN_HUM;
    }

    if state.fresh_organic_matter[layer][column] <= 0.0
        && state.humus_organic_matter[layer][column] <= 0.0
    {
        return;
    }

    let mut cn_ratio_effect = 1.0;
    let total_soil_n = state.fresh_organic_nitrogen[layer][column]
        + state.vol_no3_n_content[layer][column]
        + state.vol_nh4_n_content[layer][column];
    if total_soil_n > 0.0 {
        let cn_ratio = state.fresh_organic_matter[layer][column] * 0.4 / total_soil_n;
        if cn_ratio >= 1000.0 {
            cn_ratio_effect = 0.0;
        } else if cn_ratio > CN_MAX {
            cn_ratio_effect = (-CPAR_CN_RF * (cn_ratio - CN_MAX) / CN_MAX).exp();
        }
    }

    let wf = soil_water_effect(state, layer, column, 0.5);
    let tfac = soil_temperature_effect(state.soil_temp_daily_avrg[layer][column] - 273.161);

    let gross_release_n;
    let immobilization_rate_n;
    if state.fresh_organic_matter[layer][column] > 0.00001 {
        let g1 = tfac * wf * cn_ratio_effect * DECAY_RATE_FRESH;
        let gross_release_dw = g1 * state.fresh_organic_matter[layer][column];
        gross_release_n = g1 * state.fresh_organic_nitrogen[layer][column];

        const CPAR_N_REQ: f64 = 0.0165;
        let mut immobilization = gross_release_dw
            * (CPAR_N_REQ
                - state.fresh_organic_nitrogen[layer][column]
                    / state.fresh_organic_matter[layer][column]);
        let rnac1 = state.vol_nh4_n_content[layer][column] + state.vol_no3_n_content[layer][column]
            - 2.0 * CPAR_MIN_NH4;
        if immobilization > rnac1 {
            immobilization = rnac1;
        }
        if immobilization < 0.0 {
            immobilization = 0.0;
        }

        state.fresh_organic_matter[layer][column] -= gross_release_dw;
        state.fresh_organic_nitrogen[layer][column] += immobilization - gross_release_n;
        immobilization_rate_n = immobilization;
    } else {
        gross_release_n = 0.0;
        immobilization_rate_n = 0.0;
    }

    let rhmin = state.humus_nitrogen[layer][column] * DECAY_RATE_HUMUS * tfac * wf;
    state.humus_nitrogen[layer][column] -= rhmin + CPAR_HUMUS_N * gross_release_n;
    state.humus_organic_matter[layer][column] -=
        CN_HUM * rhmin / 0.4 + CPAR_HUMUS_N * CN_FRESH * gross_release_n / 0.4;

    let net_n_released = (1.0 - CPAR_HUMUS_N) * gross_release_n + rhmin - immobilization_rate_n;
    if net_n_released > 0.0 {
        state.vol_nh4_n_content[layer][column] += net_n_released;
        *state.mineralized_organic_n += net_n_released * state.dl[layer] * state.wk[column];
    } else {
        let mut nnom1 = 0.0;
        if state.vol_nh4_n_content[layer][column] > CPAR_MIN_NH4 {
            let addvnc =
                if net_n_released.abs() < (state.vol_nh4_n_content[layer][column] - CPAR_MIN_NH4) {
                    -net_n_released
                } else {
                    state.vol_nh4_n_content[layer][column] - CPAR_MIN_NH4
                };
            state.vol_nh4_n_content[layer][column] -= addvnc;
            *state.mineralized_organic_n -= addvnc * state.dl[layer] * state.wk[column];
            state.fresh_organic_nitrogen[layer][column] += addvnc;
            nnom1 = net_n_released + addvnc;
        }

        if nnom1 < 0.0 && state.vol_no3_n_content[layer][column] > CPAR_MIN_NH4 {
            let addvnc = if nnom1.abs() < (state.vol_no3_n_content[layer][column] - CPAR_MIN_NH4) {
                -nnom1
            } else {
                state.vol_no3_n_content[layer][column] - CPAR_MIN_NH4
            };
            state.vol_no3_n_content[layer][column] -= addvnc;
            state.fresh_organic_nitrogen[layer][column] += addvnc;
            *state.mineralized_organic_n -= addvnc * state.dl[layer] * state.wk[column];
        }
    }
}

fn soil_temperature_effect(tt: f64) -> f64 {
    let mut tfm = 0.010645 * (0.12979 * tt).exp();
    tfm = tfm.clamp(0.0, 2.0);
    tfm
}

fn nitrification(
    state: &mut SoilNitrogenState<'_>,
    layer: usize,
    column: usize,
    depth_of_layer: f64,
) {
    const CPAR_DEPTH: f64 = 0.45;
    const CPAR_NIT1: f64 = 24.635;
    const CPAR_NIT2: f64 = 8227.0;
    const CPAR_SANC: f64 = 204.0;

    let sanc = if state.vol_nh4_n_content[layer][column] < 0.1 {
        1.0 - (-CPAR_SANC * state.vol_nh4_n_content[layer][column]).exp()
    } else {
        1.0
    };

    let con1 = (CPAR_NIT1 - CPAR_NIT2 / state.soil_temp_daily_avrg[layer][column]).exp();
    let mut ratenit = 1.0 - (-con1).exp();
    let mut tff = (depth_of_layer - 30.0) / 30.0;
    if tff < 0.0 {
        tff = 0.0;
    }

    ratenit = ratenit * sanc * soil_water_effect(state, layer, column, 1.0) * CPAR_DEPTH.powf(tff);
    ratenit = ratenit.clamp(0.0, 0.10);

    let dnit = ratenit * state.vol_nh4_n_content[layer][column];
    state.vol_nh4_n_content[layer][column] -= dnit;
    state.vol_no3_n_content[layer][column] += dnit;
}

fn denitrification(state: &mut SoilNitrogenState<'_>, layer: usize, column: usize) {
    const CPAR_01: f64 = 24.5;
    const CPAR_02: f64 = 3.1;
    const CPAR_DENIT: f64 = 0.00006;
    const CPAR_FT: f64 = 0.046;
    const CPAR_HUM: f64 = 0.58;
    const VNO3_MIN: f64 = 0.00025;

    let soilc = CPAR_HUM * state.humus_organic_matter[layer][column];
    let cw = CPAR_01 + CPAR_02 * soilc;

    let mut fw = (state.vol_water_content[layer][column] - state.field_capacity[layer])
        / (state.thts[layer] - state.field_capacity[layer]);
    if fw < 0.0 {
        fw = 0.0;
    }

    let mut ft = 0.1 * (CPAR_FT * (state.soil_temp_daily_avrg[layer][column] - 273.161)).exp();
    if ft > 1.0 {
        ft = 1.0;
    }

    let mut dnrate = CPAR_DENIT * cw * state.vol_no3_n_content[layer][column] * fw * ft;
    if dnrate > state.vol_no3_n_content[layer][column] - VNO3_MIN {
        dnrate = state.vol_no3_n_content[layer][column] - VNO3_MIN;
    }
    if dnrate < 0.0 {
        dnrate = 0.0;
    }

    state.vol_no3_n_content[layer][column] -= dnrate;
    *state.soil_nitrogen_loss += dnrate * state.dl[layer] * state.wk[column];
}

pub fn soil_nitrogen_bal() {
    unsafe {
        let balsn = SoilNitrogenAtStart + CumFertilizerN + MineralizedOrganicN
            - CumNitrogenUptake
            - TotalSoilNitrogen
            - SoilNitrogenLoss;
        let _ = balsn;
    }
}

pub fn soil_nitrogen_average() {
    unsafe {
        let mut avno30 = 0.0;
        let mut avno60 = 0.0;
        let mut avno90 = 0.0;
        let mut avno120 = 0.0;
        let mut avnh30 = 0.0;
        let mut avnh60 = 0.0;
        let mut avnh90 = 0.0;
        let mut avnh120 = 0.0;

        for column in 0..nk as usize {
            for layer in 0..8 {
                avno30 += VolNo3NContent[layer][column] * dl[layer];
                avnh30 += VolNh4NContent[layer][column] * dl[layer];
            }
            for layer in 8..14 {
                avno60 += VolNo3NContent[layer][column] * dl[layer];
                avnh60 += VolNh4NContent[layer][column] * dl[layer];
            }
            for layer in 14..20 {
                avno90 += VolNo3NContent[layer][column] * dl[layer];
                avnh90 += VolNh4NContent[layer][column] * dl[layer];
            }
            for layer in 20..26 {
                avno120 += VolNo3NContent[layer][column] * dl[layer];
                avnh120 += VolNh4NContent[layer][column] * dl[layer];
            }
        }

        avno30 = 1000.0 * avno30 / (30.0 * nk as f64);
        avnh30 = 1000.0 * avnh30 / (30.0 * nk as f64);
        avno60 = 1000.0 * avno60 / (30.0 * nk as f64);
        avnh60 = 1000.0 * avnh60 / (30.0 * nk as f64);
        avno90 = 1000.0 * avno90 / (30.0 * nk as f64);
        avnh90 = 1000.0 * avnh90 / (30.0 * nk as f64);
        avno120 = 1000.0 * avno120 / (30.0 * nk as f64);
        avnh120 = 1000.0 * avnh120 / (30.0 * nk as f64);

        let _ = (
            avno30, avnh30, avno60, avnh60, avno90, avnh90, avno120, avnh120,
        );
    }
}
