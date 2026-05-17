use crate::input_functions::form;
use crate::LegacyGlobalState;
use std::os::raw::c_int;
use std::sync::{LazyLock, RwLock};

#[derive(Debug, Clone, Copy)]
struct EmergenceScratch {
    delay_of_emergence: f64,
    hypocotyl_length: f64,
    seed_moisture: f64,
    seed_layer: i32,
}

impl Default for EmergenceScratch {
    fn default() -> Self {
        Self {
            delay_of_emergence: 0.0,
            hypocotyl_length: 0.3,
            seed_moisture: 8.0,
            seed_layer: 0,
        }
    }
}

static EMERGENCE_SCRATCH: LazyLock<RwLock<EmergenceScratch>> =
    LazyLock::new(|| RwLock::new(EmergenceScratch::default()));

fn read_emergence_scratch() -> EmergenceScratch {
    *EMERGENCE_SCRATCH
        .read()
        .expect("emergence scratch state lock should not be poisoned")
}

fn write_emergence_scratch(scratch: EmergenceScratch) {
    *EMERGENCE_SCRATCH
        .write()
        .expect("emergence scratch state lock should not be poisoned") = scratch;
}

pub(crate) fn reset_scratch_state() {
    write_emergence_scratch(EmergenceScratch::default());
}

fn mark_simulation_end() -> i32 {
    let mut legacy = LegacyGlobalState::from_globals();
    legacy.b_end = true;
    let daynum = legacy.daynum;
    legacy.write_to_globals();
    daynum
}

pub(crate) fn sensible_heat_transfer(tsf: f64, tenviron: f64, height: f64, wndcanp: f64) -> f64 {
    const GRAV: f64 = 980.0;
    const S40: f64 = 0.13;
    const S42: f64 = 0.63;
    const STMIN: f64 = 5.0;
    const VONKAR: f64 = 0.40;
    const ZALIT1: f64 = 0.0962;

    let mut wind = wndcanp;
    if wind < 100.0 {
        wind = 100.0;
    }

    let mut z0 = S40 * height;
    if z0 < 1.0 {
        z0 = 1.0;
    }

    let gtop = ((200.0 - S42 * height) / z0).ln();
    let dt = tsf - tenviron;

    let mut thstar;
    let mut ustar;
    if dt >= 0.0 {
        ustar = 1.873 + 0.570172 * dt + 0.07438568 * wind;
        thstar = -0.05573 * dt + 1.987 / wind - 6.657 * dt / wind;
    } else {
        ustar = -4.4017 + 1.067 * dt + 0.25957 * wind - 0.001683 * dt * wind;
        if ustar < 5.0 {
            ustar = 5.0;
        }
        thstar = -0.0096 - 0.1149 * dt + 0.0000377 * wind + 0.0002367 * dt * wind;
        if thstar < 0.03 {
            thstar = 0.03;
        }
    }

    let mut tbot1 = tsf;
    let mut mtest: i64 = 0;
    let mut g1 = 0.0;

    loop {
        let mut previous_thstar = 0.0;
        let mut previous_ustar = 0.0;
        let mut previous_ug1 = 0.0;
        if mtest > 0 {
            tbot1 = tsf + ZALIT1 * thstar * (ustar * z0 / 15.0).powf(0.45) / VONKAR;
            previous_ustar = ustar;
            previous_thstar = thstar;
            if g1 != 0.0 {
                previous_ug1 = ustar / g1;
            }
        }

        let mut zl;
        if thstar.abs() < 1e-30 {
            zl = 0.0;
        } else {
            let mean_temp = (tenviron + tbot1) * 0.5;
            let lstar = (mean_temp * ustar * ustar) / (VONKAR * GRAV * thstar);
            zl = (200.0 - S42 * height) / lstar;
            if zl < -5.0 {
                zl = -5.0;
            }
            if zl > 0.5 {
                zl = 0.5;
            }
        }

        let (g1u, mut g2) = if zl > 0.0 {
            let g1u = -4.7 * zl;
            let mut g2 = -6.35135 * zl;
            if g2 < -1.0 {
                g2 = -1.0;
            }
            (g1u, g2)
        } else {
            let tmp1 = (1.0 - 15.0 * zl).powf(0.25);
            let g1u = 2.0 * ((1.0 + tmp1) / 2.0).ln() + ((1.0 + tmp1 * tmp1) / 2.0).ln()
                - 2.0 * (tmp1 + 1.5708).atan();
            let g2 = 2.0 * ((1.0 + (1.0 - 9.0 * zl).sqrt()) / 2.0).ln();
            (g1u, g2)
        };

        if g2 > gtop {
            g2 = gtop;
        }

        ustar = VONKAR * wind / (gtop - g1u);
        if ustar < STMIN {
            ustar = STMIN;
        }
        g1 = 0.74 * (gtop - g2) + ZALIT1 * (ustar * z0 / 0.15).powf(0.45);
        thstar = -dt * VONKAR / g1;

        if mtest > 30 {
            thstar = (thstar + previous_thstar) * 0.5;
            ustar = (ustar + previous_ustar) * 0.5;
        }

        mtest += 1;
        if mtest > 100 {
            eprintln!(" Infinite loop in SensibleHeatTransfer(). Abnormal stop!! ");
            eprintln!(" tenviron = {:10.3}", tenviron);
            eprintln!(" tsf      = {:10.3}", tsf);
            eprintln!(" PlantHeight = {:10.3}", height);
            eprintln!(" u = {:10.3}", wind);
            let _ = mark_simulation_end();
            return 0.0;
        }

        let ug1 = ustar / g1;
        let ug1res = if previous_ug1.abs() <= 1e-30 {
            ug1.abs()
        } else {
            ((previous_ug1 - ug1) / previous_ug1).abs()
        };

        if (ug1 - previous_ug1).abs() <= 0.05 || ug1res <= 0.01 {
            return ustar * VONKAR / g1;
        }
    }
}

pub(crate) fn soil_surface_balance(
    ihr: c_int,
    k: c_int,
    ess: f64,
    rlzero: f64,
    rss: f64,
    sf: f64,
    hsg: f64,
    so: &mut f64,
    so2: &mut f64,
    so3: &mut f64,
    thet: f64,
    tm: f64,
    tv: f64,
    mulch_tran_lw: f64,
    dl0: f64,
    dl1: f64,
    dl2: f64,
    vol_water0: f64,
    vol_water1: f64,
    vol_water2: f64,
    pore_space: &[f64],
    field_capacity: &[f64],
    sand_volume_fraction: &[f64],
    clay_volume_fraction: &[f64],
    marginal_water_content: &[f64],
    heat_cond_dry_soil: &[f64],
    dsand_value: f64,
    dclay_value: f64,
) {
    const EF: f64 = 0.95;
    const EG: f64 = 0.95;
    const STEFA1: f64 = 1.38e-12;

    let mut so_value = *so;
    let mut so2_value = *so2;
    let mut so3_value = *so3;

    let mut rls1 = if sf >= 0.05 {
        (1.0 - sf) * EG * rlzero + sf * EG * EF * STEFA1 * tv.powi(4)
    } else {
        EG * rlzero
    };
    if tm > 0.0 {
        rls1 = rls1 * mulch_tran_lw + EG * (1.0 - mulch_tran_lw) * STEFA1 * tm.powi(4);
    }

    let rls4 = EG * STEFA1;
    let mut previous_adjustment = 0.0;
    let mut mon = 0;
    let mut previous_so = so_value;

    while mon < 50 {
        let hlat = (75.5255 - 0.05752 * so_value) * ess;
        let dhlat = -0.05752 * ess;

        let rosoil1 = thermal_cond_soil(
            vol_water0,
            so_value,
            1,
            pore_space,
            field_capacity,
            sand_volume_fraction,
            clay_volume_fraction,
            marginal_water_content,
            heat_cond_dry_soil,
            dsand_value,
            dclay_value,
        );
        let rosoil2 = thermal_cond_soil(
            vol_water1,
            so2_value,
            2,
            pore_space,
            field_capacity,
            sand_volume_fraction,
            clay_volume_fraction,
            marginal_water_content,
            heat_cond_dry_soil,
            dsand_value,
            dclay_value,
        );
        let rosoil3 = thermal_cond_soil(
            vol_water2,
            so3_value,
            3,
            pore_space,
            field_capacity,
            sand_volume_fraction,
            clay_volume_fraction,
            marginal_water_content,
            heat_cond_dry_soil,
            dsand_value,
            dclay_value,
        );
        let rosoil = (rosoil1 * dl0 + rosoil2 * dl1 + rosoil3 * dl2)
            / (dl0 + dl1 + dl2)
            / (0.5 * dl0 + dl1 + 0.5 * dl2);
        let bbsoil = rosoil * (so_value - so3_value);
        let emtlw = rls4 * so_value.powi(4);

        let (senheat, dsenheat) = if tm > 0.0 {
            (hsg * (so_value - tm), hsg)
        } else {
            let tafk = (1.0 - sf) * thet + sf * (0.1 * so_value + 0.3 * thet + 0.6 * tv);
            (hsg * (so_value - tafk), hsg * (1.0 - sf * 0.1))
        };

        let bb = emtlw - rls1 + bbsoil + hlat - rss + senheat;
        if bb.abs() < 10e-6 {
            *so = so_value;
            *so2 = so2_value;
            *so3 = so3_value;
            return;
        }

        let demtlw = 4.0 * rls4 * so_value.powi(3);
        let rosoil1p = thermal_cond_soil(
            vol_water0,
            so_value + 0.001,
            1,
            pore_space,
            field_capacity,
            sand_volume_fraction,
            clay_volume_fraction,
            marginal_water_content,
            heat_cond_dry_soil,
            dsand_value,
            dclay_value,
        );
        let rosoilp = (rosoil1p * dl0 + rosoil2 * dl1 + rosoil3 * dl2)
            / (dl0 + dl1 + dl2)
            / (0.5 * dl0 + dl1 + 0.5 * dl2);
        let drosoil = (rosoilp - rosoil) / 0.001;
        let dbbsoil = rosoil + drosoil * (so_value - so3_value);

        let bbp = demtlw + dbbsoil + dhlat + dsenheat;
        let mut bbadjust = bb / bbp;
        if bbadjust.abs() < 0.002 {
            *so = so_value;
            *so2 = so2_value;
            *so3 = so3_value;
            return;
        }

        if mon >= 2
            && (bbadjust + previous_adjustment).abs() < (bbadjust - previous_adjustment).abs()
        {
            bbadjust = (bbadjust + previous_adjustment) / 2.0;
            so_value = (so_value + previous_so) / 2.0;
        }

        bbadjust = bbadjust.clamp(-10.0, 10.0);
        so_value -= bbadjust;
        so2_value += (so_value - previous_so) / 2.0;
        so3_value += (so_value - previous_so) / 3.0;
        previous_so = so_value;
        previous_adjustment = bbadjust;
        mon += 1;
    }

    *so = so_value;
    *so2 = so2_value;
    *so3 = so3_value;

    eprintln!(" Infinite loop in SoilSurfaceBalance(). Abnormal stop!! ");
    let daynum = mark_simulation_end();
    eprintln!("Daynum, ihr, k = {:3} {:3} {:3}", daynum, ihr, k);
    eprintln!(" so      = {:10.0}", so_value);
    eprintln!(" so2 = {:10.0}", so2_value);
    eprintln!(" so3 = {:10.0}", so3_value);
}

pub(crate) fn canopy_balance(
    ihr: c_int,
    k: c_int,
    etp1: f64,
    rlzero: f64,
    rsv: f64,
    c2: f64,
    sf: f64,
    so: f64,
    thet: f64,
    tm: f64,
    tv: &mut f64,
    mulch_tran_lw: f64,
) {
    const EF: f64 = 0.95;
    const EG: f64 = 0.95;
    const STEFA1: f64 = 1.38e-12;
    let rlv1 = if tm > 0.0 {
        sf * EF * rlzero
            + sf * EF * EG * STEFA1 * mulch_tran_lw * so.powi(4)
            + sf * EF * (1.0 - mulch_tran_lw) * STEFA1 * tm.powi(4)
    } else {
        sf * EF * rlzero + sf * EF * EG * STEFA1 * so.powi(4)
    };

    let corr = 1.0 + EG / (EF + EG - EF * EG);
    let rlv4 = STEFA1 * sf * EF * corr;
    let dsenfheat = c2 * (1.0 - 0.6 * sf);

    let mut mot = 0;
    let mut ccadx = 0.0;

    while mot < 50 {
        let hvlat = (75.5255 - 0.05752 * *tv) * etp1;
        let dhvlat = -0.05752 * etp1;
        let cclwe = rlv4 * (*tv).powi(4);
        let tvex = *tv;

        let tafk = if tm > 0.0 {
            (1.0 - sf) * thet + sf * (0.1 * tm + 0.3 * thet + 0.6 * *tv)
        } else {
            (1.0 - sf) * thet + sf * (0.1 * so + 0.3 * thet + 0.6 * *tv)
        };

        let senfheat = c2 * (*tv - tafk);
        let cc = cclwe - rlv1 + hvlat - rsv + senfheat;
        if cc.abs() < 10e-6 {
            *tv = 0.5 * (tvex + *tv);
            return;
        }

        let dcclwe = 4.0 * rlv4 * (*tv).powi(3);
        let ccp = dcclwe + dhvlat + dsenfheat;
        let mut ccadjust = cc / ccp;
        if ccadjust.abs() < 0.002 {
            return;
        }

        if mot >= 2 && (ccadjust - ccadx).abs() > (ccadjust + ccadx).abs() {
            ccadjust = 0.5 * (ccadjust + ccadx);
            *tv = 0.5 * (*tv + tvex);
        }

        ccadjust = ccadjust.clamp(-10.0, 10.0);
        *tv -= ccadjust;
        ccadx = ccadjust;
        mot += 1;
    }

    eprintln!(" Infinite loop in CanopyBalance(). Abnormal stop!! ");
    let daynum = mark_simulation_end();
    eprintln!(" Daynum, ihr, k = {:3} {:3} {:3}", daynum, ihr, k);
    eprintln!(" so      = {:10.3}", so);
    eprintln!(" tv = {:10.3}", *tv);
}

pub(crate) fn mulch_surface_balance(
    ihr: c_int,
    k: c_int,
    rlsp: f64,
    rls5: f64,
    rsm: f64,
    sf: f64,
    hsgp: f64,
    hsgm: f64,
    so: f64,
    thet: f64,
    tm: &mut f64,
    tv: f64,
) {
    let dsenheat = if sf <= 0.05 {
        hsgp * (1.0 - sf * 0.1) + hsgm
    } else {
        hsgp + hsgm
    };

    let mut mop = 0;
    let mut tmex = *tm;
    let mut tmbalex = 0.0;
    loop {
        let emtlw = rls5 * (*tm).powi(4);
        let tafk = if sf >= 0.05 {
            (1.0 - sf) * thet + sf * (0.1 * *tm + 0.3 * thet + 0.6 * tv)
        } else {
            thet
        };

        let senheat = hsgp * (*tm - tafk) + hsgm * (*tm - so);
        let tmbal = emtlw - rlsp - rsm + senheat;
        if tmbal.abs() < 10e-6 {
            return;
        }

        let demtlw = 4.0 * rls5 * (*tm).powi(3);
        let tmbalp = demtlw + dsenheat;
        let mut tmbaladjust = tmbal / tmbalp;
        if tmbaladjust.abs() < 0.002 {
            return;
        }

        if mop >= 2 && (tmbaladjust + tmbalex).abs() < (tmbaladjust - tmbalex).abs() {
            tmbaladjust = (tmbaladjust + tmbalex) / 2.0;
            *tm = (*tm + tmex) / 2.0;
        }

        tmbaladjust = tmbaladjust.clamp(-10.0, 10.0);
        *tm -= tmbaladjust;
        tmex = *tm;
        tmbalex = tmbaladjust;
        mop += 1;

        if mop >= 50 {
            break;
        }
    }

    eprintln!(" Infinite loop in MulchSurfaceBalance(). Abnormal stop!! ");
    let daynum = mark_simulation_end();
    eprintln!(" Daynum, ihr, k = {:3} {:3} {:3}", daynum, ihr, k);
    eprintln!(" so      = {:10.3}", so);
    eprintln!(" tv = {:10.3}", tv);
}

pub(crate) fn thermal_cond_soil(
    q0: f64,
    t0: f64,
    l0: c_int,
    pore_space_profile: &[f64],
    field_capacity_profile: &[f64],
    sand_volume_fraction: &[f64],
    clay_volume_fraction: &[f64],
    marginal_water_content_profile: &[f64],
    heat_cond_dry_soil_profile: &[f64],
    dsand_value: f64,
    dclay_value: f64,
) -> f64 {
    const BCLAY: f64 = 7.0;
    const BSAND: f64 = 20.0;
    const CKA: f64 = 0.0615;
    const CKW: f64 = 1.45;

    let layer = l0 as usize;
    let tcel = t0 - 273.161;
    let (
        pore_space,
        field_capacity,
        sand_fraction,
        clay_fraction,
        marginal_water_content,
        heat_cond_dry,
    ) = (
        pore_space_profile[layer],
        field_capacity_profile[layer],
        sand_volume_fraction[layer],
        clay_volume_fraction[layer],
        marginal_water_content_profile[layer],
        heat_cond_dry_soil_profile[layer],
    );

    let bb = if tcel <= 36.0 {
        0.06188
    } else if tcel <= 40.0 {
        0.06188 + (tcel - 36.0) * (0.05790 - 0.06188) / (40.0 - 36.0)
    } else {
        0.05790
    };

    let cpn = CKA + 0.05 * (bb * tcel).exp();

    let mut xair = pore_space - q0;
    if xair < 0.0 {
        xair = 0.0;
    }

    let hcond = if q0 >= field_capacity {
        let ga = 0.333 - 0.061 * xair / pore_space;
        let dair = form(cpn, CKW, ga);
        (q0 * CKW
            + dsand_value * BSAND * sand_fraction
            + dclay_value * BCLAY * clay_fraction
            + dair * cpn * xair)
            / (q0 + dsand_value * sand_fraction + dclay_value * clay_fraction + dair * xair)
    } else {
        let mut qq = q0;
        if qq < marginal_water_content {
            qq = marginal_water_content;
        }

        let ckn = CKA + (cpn - CKA) * qq / field_capacity;
        let ga = 0.041
            + 0.244 * (qq - marginal_water_content) / (field_capacity - marginal_water_content);
        let dair = form(ckn, CKW, ga);
        let mut hcond = (qq * CKW
            + dsand_value * BSAND * sand_fraction
            + dclay_value * BCLAY * clay_fraction
            + dair * ckn * xair)
            / (qq + dsand_value * sand_fraction + dclay_value * clay_fraction + dair * xair);

        if qq <= marginal_water_content {
            hcond = (hcond - heat_cond_dry) * q0 / marginal_water_content + heat_cond_dry;
        }
        hcond
    };

    hcond / 1000.0
}

pub(crate) fn predict_emergence(hour: c_int) {
    const DPL: f64 = 5.0;

    let mut legacy = LegacyGlobalState::from_globals();
    let mut scratch = read_emergence_scratch();
    let daynum = legacy.daynum;
    let dayplant = legacy.day_plant;
    let num_layers = legacy.nl as usize;
    let plant_col = legacy.plant_row_column as usize;
    let mut delay = scratch.delay_of_emergence;
    let mut hypocotyl = scratch.hypocotyl_length;
    let mut seed_moisture = scratch.seed_moisture;
    let mut seed_layer = scratch.seed_layer;

    if daynum == dayplant && hour == 0 {
        delay = 0.0;
        hypocotyl = 0.3;
        seed_moisture = 8.0;

        let mut sumdl = 0.0;
        for layer in 0..num_layers {
            sumdl += legacy.dl[layer];
            if sumdl >= DPL {
                seed_layer = layer as i32;
                break;
            }
        }
    }

    let seed_layer_index = seed_layer.max(0) as usize;
    let psi = legacy.soil_psi[[seed_layer_index, plant_col]];
    let mut te = legacy.soil_temp[[seed_layer_index, plant_col]] - 273.161;
    if te < 10.0 {
        te = 10.0;
    }

    let mut emerged = false;
    if seed_moisture <= 80.0 {
        let mut xkl = 0.0338 + 0.0000855 * te * te - 0.003479 * psi;
        if xkl < 0.0 {
            xkl = 0.0;
        }

        let dw = xkl * (80.0 - seed_moisture);
        let mut delw = if te < 21.2 {
            -0.1133 + 0.000705 * te * te - 0.001348 * psi + 0.001177 * psi * psi
        } else if te < 26.66 {
            -0.3584 + 0.001383 * te * te - 0.03509 * psi + 0.003507 * psi * psi
        } else if te < 32.3 {
            -0.6955 + 0.001962 * te * te - 0.08335 * psi + 0.007627 * psi * psi
                - 0.006411 * psi * te
        } else {
            3.3929 - 0.00197 * te * te - 0.36935 * psi + 0.00865 * psi * psi + 0.007306 * psi * te
        };

        if delw < 0.01 {
            delw = 0.01;
        }

        if dw > delw {
            seed_moisture += dw;
        } else {
            seed_moisture = 100.0;
        }
    } else {
        let xt = if te > 39.9 {
            0.0
        } else {
            0.0853 - 0.0057 * (te - 34.44) * (te - 34.44) / (41.9 - te)
        };

        if xt < 0.0 && te < 14.0 {
            delay += xt / 2.0;
        } else {
            let mut xt_adjusted = xt;
            let mut can_grow = true;
            if delay < 0.0 {
                if delay + xt_adjusted < 0.0 {
                    delay += xt_adjusted;
                    can_grow = false;
                } else {
                    xt_adjusted += delay;
                    delay = 0.0;
                }
            }

            if can_grow {
                let de = 0.0567 * xt_adjusted * hypocotyl * (10.0 - hypocotyl);
                hypocotyl += de;
                if hypocotyl > DPL {
                    emerged = true;
                }
            }
        }
    }

    scratch.delay_of_emergence = delay;
    scratch.hypocotyl_length = hypocotyl;
    scratch.seed_moisture = seed_moisture;
    scratch.seed_layer = seed_layer;
    if emerged {
        legacy.isw = 2;
        legacy.day_emerge = daynum;
        legacy.kday = 1;
    }
    legacy.write_to_globals();
    write_emergence_scratch(scratch);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn soil_profiles() -> (
        [f64; 40],
        [f64; 40],
        [f64; 40],
        [f64; 40],
        [f64; 40],
        [f64; 40],
    ) {
        let mut pore_space = [0.45; 40];
        let mut field_capacity = [0.30; 40];
        let mut sand_fraction = [0.35; 40];
        let mut clay_fraction = [0.25; 40];
        let mut marginal_water = [0.05; 40];
        let mut dry_cond = [0.08; 40];
        // Keep layer 0/1/2 explicitly initialized to meaningful values used in tests.
        for i in 0..3 {
            pore_space[i] = 0.45;
            field_capacity[i] = 0.30;
            sand_fraction[i] = 0.35;
            clay_fraction[i] = 0.25;
            marginal_water[i] = 0.05;
            dry_cond[i] = 0.08;
        }
        (
            pore_space,
            field_capacity,
            sand_fraction,
            clay_fraction,
            marginal_water,
            dry_cond,
        )
    }

    #[test]
    fn thermal_cond_soil_is_positive_and_higher_when_wetter() {
        let (pore_space, field_capacity, sand_fraction, clay_fraction, marginal_water, dry_cond) =
            soil_profiles();

        let dry = thermal_cond_soil(
            0.06,
            293.161,
            1,
            &pore_space,
            &field_capacity,
            &sand_fraction,
            &clay_fraction,
            &marginal_water,
            &dry_cond,
            0.2,
            0.2,
        );
        let wet = thermal_cond_soil(
            0.25,
            293.161,
            1,
            &pore_space,
            &field_capacity,
            &sand_fraction,
            &clay_fraction,
            &marginal_water,
            &dry_cond,
            0.2,
            0.2,
        );

        assert!(dry.is_finite() && wet.is_finite());
        assert!(dry > 0.0);
        assert!(wet > dry, "expected wetter soil conductivity to be higher");
    }

    #[test]
    fn sensible_heat_transfer_returns_finite_positive_value() {
        let coeff = sensible_heat_transfer(305.0, 300.0, 50.0, 180.0);
        assert!(coeff.is_finite());
        assert!(coeff > 0.0);
    }
}
