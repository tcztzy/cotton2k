use chrono::prelude::*;
use chrono::Duration;
use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

#[cfg(test)]
mod tests;

mod root;
mod stem;

#[no_mangle]
pub extern "C" fn LeapYear(year: u64) -> u64 {
    let is_divisible = |n| year % n == 0;

    (is_divisible(4) && (!is_divisible(100) || is_divisible(400))) as u64
}

#[no_mangle]
pub extern "C" fn DateToDoy(date_str: *const c_char, year_start: i32) -> i64 {
    let date_string = unsafe { CStr::from_ptr(date_str) };
    let date = date_string.to_str().unwrap().trim();
    match NaiveDate::parse_from_str(date, "%d-%b-%Y") {
        Ok(nd) => {
            if nd.year() == year_start + 1 {
                (nd.ordinal() + NaiveDate::from_ymd(nd.year(), 1, 1).pred().ordinal()) as i64
            } else {
                nd.ordinal() as i64
            }
        }
        Err(..) => 0,
    }
}

#[no_mangle]
pub extern "C" fn DoyToDate(doy: i32, year_start: i32) -> *const u8 {
    if doy > 0 {
        (NaiveDate::from_ymd(year_start, 1, 1).pred() + Duration::days(doy.into()))
            .format("%d-%b-%Y\0")
            .to_string()
            .to_uppercase()
            .as_ptr()
    } else {
        (" ".repeat(11) + "\0").as_ptr()
    }
}

/// This function sorts an array of values by its value (larger is first)
/// together with three indexes associated with each value.
#[no_mangle]
pub extern "C" fn SortArray(size: usize, data: *mut f64, ik: *mut i32, il: *mut i32, im: *mut i32) {
    let _data = unsafe { slice::from_raw_parts_mut(data, size) };
    let _ik = unsafe { slice::from_raw_parts_mut(ik, size) };
    let _il = unsafe { slice::from_raw_parts_mut(il, size) };
    let _im = unsafe { slice::from_raw_parts_mut(im, size) };
    let mut x = [(0f64, 0i32, 0i32, 0i32)].repeat(size);
    for i in 0..size {
        x[i] = (_data[i], _ik[i], _il[i], _im[i]);
    }
    x.sort_by(|a, b| (b.0).partial_cmp(&a.0).unwrap());
    unsafe {
        for i in 0..size {
            *(data.offset(i as isize)) = x[i].0;
            *(ik.offset(i as isize)) = x[i].1;
            *(il.offset(i as isize)) = x[i].2;
            *(im.offset(i as isize)) = x[i].3;
        }
    }
}

///     This function computes and returns the resistance of leaves of cotton
/// plants to transpiration. It is assumed to be a function of leaf age.
/// It is called from LeafWaterPotential.
///     The input argument (age) is leaf age in physiological days.
#[no_mangle]
pub extern "C" fn LeafResistance(age: f64) -> f64 {
    // The following constant parameters are used:
    let afac: f64 = 160.; // factor used for computing leaf resistance.
    let agehi: f64 = 94.; // higher limit for leaf age.
    let agelo: f64 = 48.; // lower limit for leaf age.
    let rlmin: f64 = 0.5; // minimum leaf resistance.

    if age <= agelo {
        rlmin
    } else if age >= agehi {
        rlmin + (agehi - agelo) * (agehi - agelo) / afac
    } else {
        let ax: f64 = 2. * agehi - agelo; // intermediate variable
        rlmin + (age - agelo) * (ax - age) / afac
    }
}

/// computes physiological age
///     This function returns the daily 'physiological age' increment,
///  based on hourly temperatures. It is called each day by SimulateThisDay.
#[no_mangle]
pub extern "C" fn PhysiologicalAge(air_temp: *const f64) -> f64 {
    let air_temperature = unsafe { slice::from_raw_parts(air_temp, 24) };
    //     The threshold value is assumed to be 12 C (p1). One physiological day is
    //  equivalent to a day with an average temperature of 26 C, and therefore the
    //  heat units are divided by 14 (p2).
    //     A linear relationship is assumed between temperature and heat unit
    //  accumulation in the range of 12 C (p1) to 33 C (p2*p3+p1). the effect of
    //  temperatures higher than 33 C is assumed to be equivalent to that of 33 C.
    //     The following constant Parameters are used in this function:
    let p1 = 12.; // threshold temperature, C
    let p2 = 14.; // temperature, C, above p1, for one physiological day.
    let p3 = 1.5; // maximum value of a physiological day.

    let mut dayfd = 0.; // the daily contribution to physiological age (return value).
    for t in air_temperature {
        let mut tfd = (t - p1) / p2; // the hourly contribution to physiological age.
        if tfd < 0. {
            tfd = 0.;
        }
        if tfd > p3 {
            tfd = p3;
        }
        dayfd += tfd;
    }
    dayfd / 24.
}

/// This function is called from PotentialRootGrowth(), TapRootGrowth() and
/// LateralRootGrowth(). It computes the effects of soil temperature on the rate
/// growth. It is essentially based on the usage of GOSSYM, but relative values
/// are computed here. The computed value returned by this function is between 0 and 1.
///
/// It is assumed that maximum root growth occurs at or above 30 C, and no root growth
/// occurs at or below 13.5 C. A quadratic response to temperature between these limits
/// is assumed.
///
/// The following argument is used:
///
/// t - Soil temperature (C), daily average.
///       
/// The parameters used are p1, p2, p3, with the following results:
/// t =      14    16    18    20    22    24    26    28    30   
/// trf =  .053  .261  .443  .600  .731  .837  .917  .971  1.00
///
#[no_mangle]
pub extern "C" fn SoilTemOnRootGrowth(t: f64) -> f64 {
    let p1 = -2.12;
    let p2 = 0.2;
    let p3 = -0.0032;

    if t >= 30. {
        1.
    } else {
        let trf = p1 + t * (p2 + p3 * t);
        if trf > 1. {
            1.
        } else if trf < 0. {
            0.
        } else {
            trf
        }
    }
}

/// effects of pix applied.
#[no_mangle]
pub extern "C" fn Pix() { // TODO
}

///  This is the temperature function for leaf growth rate. It is called by function
///  PotentialLeafGrowth(). It is based on the original code of GOSSYM, and the parameters
///  are the same. The argument t is air temperature (C).
///
///  ra is divided by the maximum value. Thus, the function returns values between `0.` and `1.`
///
///  This will result :
///
///     maximum value of TemperatureOnLeafGrowthRate = 1.00  at  t = 29.86747
///         for t = 24    TemperatureOnLeafGrowthRate = 0.766
///         for t = 27    TemperatureOnLeafGrowthRate = 0.953
///         for t = 30    TemperatureOnLeafGrowthRate = 0.999
///         for t = 36    TemperatureOnLeafGrowthRate = 0.737
///         for t = 42    TemperatureOnLeafGrowthRate = 0.0
///      and for t <= 24 :
///         for t = 12    TemperatureOnLeafGrowthRate = 0.0
///         for t = 16    TemperatureOnLeafGrowthRate = 0.269
///         for t = 20    TemperatureOnLeafGrowthRate = 0.549
///         for t = 24    TemperatureOnLeafGrowthRate = 0.768
#[no_mangle]
pub extern "C" fn TemperatureOnLeafGrowthRate(t: f64) -> f64 {
    // Constant parameters used:
    let par: [f64; 8] = [
        24.,
        -1.14277,
        0.0910026,
        0.00152344,
        -0.317136,
        0.0300712,
        0.000416356,
        0.2162044,
    ];
    // intermediate value for computing temperatureOnLeafGrowthRate.
    let ra = if t > par[0] {
        par[1] + t * (par[2] - t * par[3])
    } else {
        par[4] + t * (par[5] - t * par[6])
    };
    if ra < 0. {
        0.
    } else {
        ra / par[7]
    }
}

/// This function computes the effect of air temperature (t) on growth
/// rate of bolls in cotton plants. It is called from PotentialFruitGrowth().
/// Some values computed by this function:
/// t (C)       tfr
/// 12          0.
/// 15          0.336
/// 20          0.751
/// 25          0.978
/// 26          1.
/// 28.5        1.024 (maximum)
/// 30          1.016
/// 35          0.866
/// 40          0.527
/// 45          0.
#[no_mangle]
pub extern "C" fn TemperatureOnFruitGrowthRate(t: f64) -> f64 {
    let p1 = -2.041;
    let p2 = 0.215;
    let p3 = 0.00377;

    let tfr = p1 + t * (p2 - p3 * t);
    if tfr < 0. {
        0.
    } else {
        tfr
    }
}

/// The function PsiOnTranspiration() computes and returns the effect of the average soil
/// matrix water potential on transpiration rate. It is called by WaterUptake().
/// The argument PsiAverage is the average soil water matrix potential, bars.
#[no_mangle]
pub extern "C" fn PsiOnTranspiration(psi_average: f64) -> f64 {
    // This is a third degree function with two parameters (a, b). It has the
    // value of 1 when PsiAverage = b - a, and the value of 0 when PsiAverage = - a.

    // The minimum value, however, is set to d, and the maximum value to c.
    let a = 20.;
    let b = 14.;
    let c = 1.00;
    let d = 0.05;
    let rfep = ((a + psi_average) / b).powi(3);
    if rfep > c {
        c
    } else if rfep < d {
        d
    } else {
        rfep
    }
}

/// This function computes soil water osmotic potential (in bars, positive value).
#[no_mangle]
pub extern "C" fn PsiOsmotic(q: f64, qsat: f64, ec: f64) -> f64
// The following arguments are used:
//   q - soil water content, cm3 cm-3.
//   qsat - saturated water content, cm3 cm-3.
//   ec - electrical conductivity of saturated extract (mmho/cm)
{
    if ec > 0f64 {
        let result = 0.36 * ec * qsat / q;
        if result > 6f64 {
            6f64
        } else {
            result
        }
    } else {
        0f64
    }
}

/// This function computes soil water content (cm3 cm-3) for
/// a given value of matrix potential, using the Van-Genuchten equation.
#[no_mangle]
pub extern "C" fn qpsi(psi: f64, qr: f64, qsat: f64, alpha: f64, beta: f64) -> f64
// The following arguments are used:
//   alpha, beta  - parameters of the van-genuchten equation.
//   psi - soil water matrix potential (bars).
//   qr - residual water content, cm3 cm-3.
//   qsat - saturated water content, cm3 cm-3.
{
    // For very high values of PSI, saturated water content is assumed.
    // For very low values of PSI, air-dry water content is assumed.
    if psi >= -0.00001 {
        qsat
    } else if psi <= -500000f64 {
        qr
    } else {
        // The soil water matric potential is transformed from bars (psi)
        // to cm in positive value (psix).
        let psix = 1000. * (psi + 0.00001).abs();
        // The following equation is used (in FORTRAN notation):
        //   QPSI = QR + (QSAT-QR) / (1 + (ALPHA*PSIX)**BETA)**(1-1/BETA)
        let gama = 1. - 1. / beta;
        let term = 1. + (alpha * psix).powf(beta); //  intermediate variable
        let swfun = qr + (qsat - qr) / term.powf(gama); //  computed water content
        if swfun < (qr + 0.0001) {
            qr + 0.0001
        } else {
            swfun
        }
    }
}

/// This function computes soil water matric potential (in bars) for a given value of soil water content, using the Van-Genuchten equation.
#[no_mangle]
pub extern "C" fn psiq(q: f64, qr: f64, qsat: f64, alpha: f64, beta: f64) -> f64
// The following arguments are used:
//   alpha, beta  - parameters of the van-genuchten equation.
//   q - soil water content, cm3 cm-3.
//   qr - residual water content, cm3 cm-3.
//   qsat - saturated water content, cm3 cm-3.
{
    // For very low values of water content (near the residual water
    // content) psiq is -500000 bars, and for saturated or higher water
    // content psiq is -0.00001 bars.
    if (q - qr) < 0.00001 {
        return -500000.;
    } else if q >= qsat {
        return -0.00001;
    }
    // The following equation is used (FORTRAN notation):
    // PSIX = (((QSAT-QR) / (Q-QR))**(1/GAMA) - 1) **(1/BETA) / ALPHA
    let gama = 1. - 1. / beta;
    let gaminv = 1. / gama;
    let term = ((qsat - qr) / (q - qr)).powf(gaminv); //  intermediate variable
    let mut psix = (term - 1.).powf(1. / beta) / alpha;
    if psix < 0.01 {
        psix = 0.01;
    }
    // psix (in cm) is converted to bars (negative value).
    psix = (0.01 - psix) * 0.001;
    if psix < -500000. {
        psix = -500000.;
    }
    if psix > -0.00001 {
        psix = -0.00001;
    }
    return psix;
}

/// This function computes soil water hydraulic conductivity
/// for a given value of soil water content, using the Van-Genuchten
/// equation. The units of the computed conductivity are the same as the given
/// saturated conductivity (SaturatedHydCond).
#[no_mangle]
pub extern "C" fn wcond(
    q: f64,
    qr: f64,
    qsat: f64,
    beta: f64,
    saturated_hyd_cond: f64,
    pore_space: f64,
) -> f64
// The following arguments are used:
//   beta  - parameter of the van-genuchten equation.
//   saturated_hyd_cond - saturated hydraulic conductivity (at qsat).
//   pore_space - pore space volume.
//   q - soil water content, cm3 cm-3.
//   qr - residual water content, cm3 cm-3.
//   qsat - saturated water content, cm3 cm-3.
{
    // For very low values of water content (near the residual water content) wcond is 0.
    if (q - qr) < 0.0001 {
        return 0.;
    }
    // Water content for saturated conductivity is minimum of PoreSpace and qsat.

    // For very high values of water content (exceeding the saturated
    // water content or pore space) conductivity is SaturatedHydCond.
    let xsat = if qsat < pore_space { qsat } else { pore_space };
    if q >= xsat {
        return saturated_hyd_cond;
    }
    // The following equation is used (in FORTRAN notation):
    //   WCOND = CONDSAT * ((Q-QR)/(XSAT-QR))**0.5
    //           * (1-(1-((Q-QR)/(XSAT-QR))**(1/GAMA))**GAMA)**2
    let gama = 1. - 1. / beta;
    let gaminv = 1. / gama;
    let sweff = (q - qr) / (xsat - qr); // intermediate variable (effective water content).
    let acoeff = (1. - sweff.powf(gaminv)).powf(gama); // intermediate variable
    let bcoeff = (1. - acoeff).powi(2); // intermediate variable
    let conductivity = sweff.powf(0.5) * bcoeff * saturated_hyd_cond;
    return conductivity;
}
