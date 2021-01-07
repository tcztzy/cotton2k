use chrono::prelude::*;
use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

#[cfg(test)]
mod tests;

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

///     This function is called from PotentialRootGrowth(), TapRootGrowth() and
///  LateralRootGrowth(). It computes the effects of soil temperature on the rate
///  growth. It is essentially based on the usage of GOSSYM, but relative values
///  are computed here. The computed value returned by this function is between 0 and 1.
///     It is assumed that maximum root growth occurs at or above 30 C, and no root growth
///  occurs at or below 13.5 C. A quadratic response to temperature between these limits
///  is assumed.
///
///     The following argument is used:
///        t - Soil temperature (C), daily average.
///       
///     The parameters used are p1, p2, p3, with the following results:
///  t =      14    16    18    20    22    24    26    28    30   
///  trf =  .053  .261  .443  .600  .731  .837  .917  .971  1.00
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
pub extern "C" fn Pix() { /* TODO */ }
