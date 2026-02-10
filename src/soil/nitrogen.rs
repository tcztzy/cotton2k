use crate::{
    dl, nk, nl, thetar, thts, wk, BulkDensity, CumFertilizerN, CumNitrogenUptake, DayStart, Daynum,
    FieldCapacity, FreshOrganicMatter, FreshOrganicNitrogen, HumusNitrogen, HumusOrganicMatter,
    MineralizedOrganicN, SoilHorizonNum, SoilNitrogenAtStart, SoilNitrogenLoss, SoilTempDailyAvrg,
    TotalSoilNitrogen, VolNh4NContent, VolNo3NContent, VolUreaNContent, VolWaterContent,
};

static mut SOIL_N_DEPTH: [f64; 40] = [0.0; 40];

pub unsafe fn soil_nitrogen() {
    if Daynum <= DayStart {
        let mut sumdl = 0.0;
        for layer in 0..nl as usize {
            sumdl += dl[layer];
            SOIL_N_DEPTH[layer] = sumdl;
        }
    }

    for layer in 0..nl as usize {
        for column in 0..nk as usize {
            if VolUreaNContent[layer][column] > 0.0 {
                urea_hydrolysis(layer, column);
            }
            mineralize_nitrogen(layer, column);
            if VolNh4NContent[layer][column] > 0.00001 {
                nitrification(layer, column, SOIL_N_DEPTH[layer]);
            }

            if VolNo3NContent[layer][column] > 0.001
                && VolWaterContent[layer][column] > FieldCapacity[layer]
                && SoilTempDailyAvrg[layer][column] >= 278.161
            {
                denitrification(layer, column);
            }
        }
    }
}

unsafe fn urea_hydrolysis(layer: usize, column: usize) {
    const AK0: f64 = 0.25;
    const CAK1: f64 = 0.3416;
    const CAK2: f64 = 0.0776;
    const STF1: f64 = 40.0;
    const STF2: f64 = 0.20;
    const SWF1: f64 = 0.20;

    let horizon = SoilHorizonNum[layer] as usize;
    let oc = 0.4 * (FreshOrganicMatter[layer][column] + HumusOrganicMatter[layer][column]) * 0.1
        / BulkDensity[horizon];

    let mut ak = CAK1 + CAK2 * oc;
    if ak < AK0 {
        ak = AK0;
    }

    let mut swf = soil_water_effect(layer, column, 0.5) + SWF1;
    swf = swf.clamp(0.0, 1.0);

    let mut stf = (SoilTempDailyAvrg[layer][column] - 273.161) / STF1 + STF2;
    stf = stf.clamp(0.0, 1.0);

    let mut hydrur = ak * swf * stf * VolUreaNContent[layer][column];
    if hydrur > VolUreaNContent[layer][column] {
        hydrur = VolUreaNContent[layer][column];
    }
    VolUreaNContent[layer][column] -= hydrur;
    VolNh4NContent[layer][column] += hydrur;
}

unsafe fn soil_water_effect(layer: usize, column: usize, xx: f64) -> f64 {
    let mut wf = if VolWaterContent[layer][column] <= FieldCapacity[layer] {
        (VolWaterContent[layer][column] - thetar[layer]) / (FieldCapacity[layer] - thetar[layer])
    } else {
        1.0 - xx * (VolWaterContent[layer][column] - FieldCapacity[layer])
            / (thts[layer] - FieldCapacity[layer])
    };

    if wf < 0.0 {
        wf = 0.0;
    }
    wf
}

unsafe fn mineralize_nitrogen(layer: usize, column: usize) {
    const CN_FRESH: f64 = 25.0;
    const CN_HUM: f64 = 10.0;
    const CN_MAX: f64 = 13.0;
    const CPAR_CN_RF: f64 = 0.693;
    const CPAR_HUMUS_N: f64 = 0.20;
    const CPAR_MIN_NH4: f64 = 0.00025;
    const DECAY_RATE_FRESH: f64 = 0.03;
    const DECAY_RATE_HUMUS: f64 = 0.000083;

    if Daynum <= DayStart {
        FreshOrganicNitrogen[layer][column] = FreshOrganicMatter[layer][column] * 0.4 / CN_FRESH;
        HumusNitrogen[layer][column] = HumusOrganicMatter[layer][column] * 0.4 / CN_HUM;
    }

    if FreshOrganicMatter[layer][column] <= 0.0 && HumusOrganicMatter[layer][column] <= 0.0 {
        return;
    }

    let mut cn_ratio_effect = 1.0;
    let total_soil_n = FreshOrganicNitrogen[layer][column]
        + VolNo3NContent[layer][column]
        + VolNh4NContent[layer][column];
    if total_soil_n > 0.0 {
        let cn_ratio = FreshOrganicMatter[layer][column] * 0.4 / total_soil_n;
        if cn_ratio >= 1000.0 {
            cn_ratio_effect = 0.0;
        } else if cn_ratio > CN_MAX {
            cn_ratio_effect = (-CPAR_CN_RF * (cn_ratio - CN_MAX) / CN_MAX).exp();
        }
    }

    let wf = soil_water_effect(layer, column, 0.5);
    let tfac = soil_temperature_effect(SoilTempDailyAvrg[layer][column] - 273.161);

    let gross_release_n;
    let immobilization_rate_n;
    if FreshOrganicMatter[layer][column] > 0.00001 {
        let g1 = tfac * wf * cn_ratio_effect * DECAY_RATE_FRESH;
        let gross_release_dw = g1 * FreshOrganicMatter[layer][column];
        gross_release_n = g1 * FreshOrganicNitrogen[layer][column];

        const CPAR_N_REQ: f64 = 0.0165;
        let mut immobilization = gross_release_dw
            * (CPAR_N_REQ
                - FreshOrganicNitrogen[layer][column] / FreshOrganicMatter[layer][column]);
        let rnac1 =
            VolNh4NContent[layer][column] + VolNo3NContent[layer][column] - 2.0 * CPAR_MIN_NH4;
        if immobilization > rnac1 {
            immobilization = rnac1;
        }
        if immobilization < 0.0 {
            immobilization = 0.0;
        }

        FreshOrganicMatter[layer][column] -= gross_release_dw;
        FreshOrganicNitrogen[layer][column] += immobilization - gross_release_n;
        immobilization_rate_n = immobilization;
    } else {
        gross_release_n = 0.0;
        immobilization_rate_n = 0.0;
    }

    let rhmin = HumusNitrogen[layer][column] * DECAY_RATE_HUMUS * tfac * wf;
    HumusNitrogen[layer][column] -= rhmin + CPAR_HUMUS_N * gross_release_n;
    HumusOrganicMatter[layer][column] -=
        CN_HUM * rhmin / 0.4 + CPAR_HUMUS_N * CN_FRESH * gross_release_n / 0.4;

    let net_n_released = (1.0 - CPAR_HUMUS_N) * gross_release_n + rhmin - immobilization_rate_n;
    if net_n_released > 0.0 {
        VolNh4NContent[layer][column] += net_n_released;
        MineralizedOrganicN += net_n_released * dl[layer] * wk[column];
    } else {
        let mut nnom1 = 0.0;
        if VolNh4NContent[layer][column] > CPAR_MIN_NH4 {
            let addvnc;
            if net_n_released.abs() < (VolNh4NContent[layer][column] - CPAR_MIN_NH4) {
                addvnc = -net_n_released;
            } else {
                addvnc = VolNh4NContent[layer][column] - CPAR_MIN_NH4;
            }
            VolNh4NContent[layer][column] -= addvnc;
            MineralizedOrganicN -= addvnc * dl[layer] * wk[column];
            FreshOrganicNitrogen[layer][column] += addvnc;
            nnom1 = net_n_released + addvnc;
        }

        if nnom1 < 0.0 && VolNo3NContent[layer][column] > CPAR_MIN_NH4 {
            let addvnc;
            if nnom1.abs() < (VolNo3NContent[layer][column] - CPAR_MIN_NH4) {
                addvnc = -nnom1;
            } else {
                addvnc = VolNo3NContent[layer][column] - CPAR_MIN_NH4;
            }
            VolNo3NContent[layer][column] -= addvnc;
            FreshOrganicNitrogen[layer][column] += addvnc;
            MineralizedOrganicN -= addvnc * dl[layer] * wk[column];
        }
    }
}

unsafe fn soil_temperature_effect(tt: f64) -> f64 {
    let mut tfm = 0.010645 * (0.12979 * tt).exp();
    tfm = tfm.clamp(0.0, 2.0);
    tfm
}

unsafe fn nitrification(layer: usize, column: usize, depth_of_layer: f64) {
    const CPAR_DEPTH: f64 = 0.45;
    const CPAR_NIT1: f64 = 24.635;
    const CPAR_NIT2: f64 = 8227.0;
    const CPAR_SANC: f64 = 204.0;

    let sanc = if VolNh4NContent[layer][column] < 0.1 {
        1.0 - (-CPAR_SANC * VolNh4NContent[layer][column]).exp()
    } else {
        1.0
    };

    let con1 = (CPAR_NIT1 - CPAR_NIT2 / SoilTempDailyAvrg[layer][column]).exp();
    let mut ratenit = 1.0 - (-con1).exp();
    let mut tff = (depth_of_layer - 30.0) / 30.0;
    if tff < 0.0 {
        tff = 0.0;
    }

    ratenit = ratenit * sanc * soil_water_effect(layer, column, 1.0) * CPAR_DEPTH.powf(tff);
    ratenit = ratenit.clamp(0.0, 0.10);

    let dnit = ratenit * VolNh4NContent[layer][column];
    VolNh4NContent[layer][column] -= dnit;
    VolNo3NContent[layer][column] += dnit;
}

unsafe fn denitrification(layer: usize, column: usize) {
    const CPAR_01: f64 = 24.5;
    const CPAR_02: f64 = 3.1;
    const CPAR_DENIT: f64 = 0.00006;
    const CPAR_FT: f64 = 0.046;
    const CPAR_HUM: f64 = 0.58;
    const VNO3_MIN: f64 = 0.00025;

    let soilc = CPAR_HUM * HumusOrganicMatter[layer][column];
    let cw = CPAR_01 + CPAR_02 * soilc;

    let mut fw = (VolWaterContent[layer][column] - FieldCapacity[layer])
        / (thts[layer] - FieldCapacity[layer]);
    if fw < 0.0 {
        fw = 0.0;
    }

    let mut ft = 0.1 * (CPAR_FT * (SoilTempDailyAvrg[layer][column] - 273.161)).exp();
    if ft > 1.0 {
        ft = 1.0;
    }

    let mut dnrate = CPAR_DENIT * cw * VolNo3NContent[layer][column] * fw * ft;
    if dnrate > VolNo3NContent[layer][column] - VNO3_MIN {
        dnrate = VolNo3NContent[layer][column] - VNO3_MIN;
    }
    if dnrate < 0.0 {
        dnrate = 0.0;
    }

    VolNo3NContent[layer][column] -= dnrate;
    SoilNitrogenLoss += dnrate * dl[layer] * wk[column];
}

pub unsafe fn soil_nitrogen_bal() {
    let balsn = SoilNitrogenAtStart + CumFertilizerN + MineralizedOrganicN
        - CumNitrogenUptake
        - TotalSoilNitrogen
        - SoilNitrogenLoss;
    let _ = balsn;
}

pub unsafe fn soil_nitrogen_average() {
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
