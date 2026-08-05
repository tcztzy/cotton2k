//! Shared scalar equations for soil water, osmotic potential, and climate lookup.
//!
//! These routines are the numerical utility boundary used by the atmosphere,
//! soil, plant, and profile modules. They transform unit-bearing model scalars
//! into model-native values and read the process-wide climate buffer when
//! requested; they do not mutate caller-owned arrays.

use crate::{
    Clim, CLIMATE_METRIC, CLIMATE_METRIC_IRRD, CLIMATE_METRIC_RAIN, CLIMATE_METRIC_TDEW,
    CLIMATE_METRIC_TMAX, CLIMATE_METRIC_TMIN, CLIMATE_METRIC_WIND,
};
use std::os::raw::c_int;

/// Converts volumetric water content to soil water potential.
pub fn psiq(q: f64, qr: f64, qsat: f64, alpha: f64, beta: f64) -> f64 {
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

/// Converts soil water potential to volumetric water content.
pub fn qpsi(psi: f64, qr: f64, qsat: f64, alpha: f64, beta: f64) -> f64 {
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

/// Computes unsaturated hydraulic conductivity from water content.
pub fn wcond(
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

/// Computes osmotic potential from water content and electrical conductivity.
pub fn PsiOsmotic(q: f64, qsat: f64, ec: f64) -> f64 {
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

/// Reads one metric for a day-of-year from the shared climate buffer.
pub fn GetFromClim(item: CLIMATE_METRIC, doy: c_int) -> f64 {
    const CLIM_LEN: usize = 400;
    let clim = Clim.read().expect("Clim lock poisoned");
    let mut i = 0usize;
    while i < CLIM_LEN {
        let record = clim[i];
        if record.nDay == doy {
            break;
        }
        i += 1;
    }
    if i >= CLIM_LEN {
        i = CLIM_LEN - 1;
    }
    if doy < clim[0].nDay {
        i = 0;
    }
    let record = clim[i];
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
