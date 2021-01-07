use chrono::prelude::*;
use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

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
/// It is called from LeafWaterPotential().
///     The input argument (age) is leaf age in physiological days.
#[no_mangle]
pub extern "C" fn LeafResistance(age: f64) -> f64 {
    // The following constant parameters are used:
    let afac: f64 = 160.;  // factor used for computing leaf resistance.
    let agehi: f64 = 94.;  // higher limit for leaf age.
    let agelo: f64 = 48.;  // lower limit for leaf age.
    let rlmin: f64 = 0.5; // minimum leaf resistance.

    if age <= agelo {
        rlmin
    }
    else if age >= agehi {
        rlmin + (agehi - agelo) * (agehi - agelo) / afac
    }
    else {
        let ax: f64 = 2. * agehi - agelo; // intermediate variable
        rlmin + (age - agelo) * (ax - age) / afac
    }
}
