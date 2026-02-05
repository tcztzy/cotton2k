use crate::{
    Clim, Climstruct, CLIMATE_METRIC, CLIMATE_METRIC_IRRD, CLIMATE_METRIC_RAIN,
    CLIMATE_METRIC_TDEW, CLIMATE_METRIC_TMAX, CLIMATE_METRIC_TMIN, CLIMATE_METRIC_WIND,
};
use std::os::raw::c_int;

#[no_mangle]
pub extern "C" fn psiq(q: f64, qr: f64, qsat: f64, alpha: f64, beta: f64) -> f64 {
    if (q - qr) < 0.00001 {
        return -500000.0;
    } else if q >= qsat {
        return -0.00001;
    }
    let gama = 1.0 - 1.0 / beta;
    let gaminv = 1.0 / gama;
    let mut term = (qsat - qr) / (q - qr);
    term = term.powf(gaminv);
    let mut psix = (term - 1.0).powf(1.0 / beta) / alpha;
    if psix < 0.01 {
        psix = 0.01;
    }
    psix = (0.01 - psix) * 0.001;
    if psix < -500000.0 {
        psix = -500000.0;
    }
    if psix > -0.00001 {
        psix = -0.00001;
    }
    psix
}

#[no_mangle]
pub extern "C" fn qpsi(psi: f64, qr: f64, qsat: f64, alpha: f64, beta: f64) -> f64 {
    if psi >= -0.00001 {
        return qsat;
    } else if psi <= -500000.0 {
        return qr;
    }
    let psix = 1000.0 * (psi + 0.00001).abs();
    let gama = 1.0 - 1.0 / beta;
    let term = 1.0 + (alpha * psix).powf(beta);
    let mut swfun = qr + (qsat - qr) / term.powf(gama);
    if swfun < (qr + 0.0001) {
        swfun = qr + 0.0001;
    }
    swfun
}

#[no_mangle]
pub extern "C" fn wcond(
    q: f64,
    qr: f64,
    qsat: f64,
    beta: f64,
    saturated_hyd_cond: f64,
    pore_space: f64,
) -> f64 {
    if (q - qr) < 0.0001 {
        return 0.0;
    }
    let xsat = qsat.min(pore_space);
    if q >= xsat {
        return saturated_hyd_cond;
    }
    let gama = 1.0 - 1.0 / beta;
    let gaminv = 1.0 / gama;
    let sweff = (q - qr) / (xsat - qr);
    let acoeff = (1.0 - sweff.powf(gaminv)).powf(gama);
    let bcoeff = (1.0 - acoeff).powi(2);
    sweff.powf(0.5) * bcoeff * saturated_hyd_cond
}

#[no_mangle]
pub extern "C" fn PsiOsmotic(q: f64, qsat: f64, ec: f64) -> f64 {
    if ec > 0.0 {
        let mut value = 0.36 * ec * qsat / q;
        if value > 6.0 {
            value = 6.0;
        }
        value
    } else {
        0.0
    }
}

#[no_mangle]
pub extern "C" fn GetFromClim(item: CLIMATE_METRIC, doy: c_int) -> f64 {
    const CLIM_LEN: usize = 400;
    unsafe {
        let clim_ptr = std::ptr::addr_of_mut!(Clim) as *mut Climstruct;
        let mut i = 0usize;
        while i < CLIM_LEN {
            let record = std::ptr::read(clim_ptr.add(i));
            if record.nDay == doy {
                break;
            }
            i += 1;
        }
        if i >= CLIM_LEN {
            i = CLIM_LEN - 1;
        }
        let first = std::ptr::read(clim_ptr);
        if doy < first.nDay {
            i = 0;
        }
        let record = std::ptr::read(clim_ptr.add(i));
        if item == CLIMATE_METRIC_TMIN {
            record.Tmin
        } else if item == CLIMATE_METRIC_TMAX {
            record.Tmax
        } else if item == CLIMATE_METRIC_IRRD {
            record.Rad
        } else if item == CLIMATE_METRIC_RAIN {
            record.Rain
        } else if item == CLIMATE_METRIC_WIND {
            record.Wind
        } else if item == CLIMATE_METRIC_TDEW {
            record.Tdew
        } else {
            -99.0
        }
    }
}
