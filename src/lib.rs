use chrono::prelude::*;
use std::ffi::CStr;
use std::os::raw::c_char;

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
