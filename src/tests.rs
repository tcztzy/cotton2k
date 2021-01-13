use super::root::*;
use super::util::*;
use std::slice;
use std::str;

fn approx_equal(a: f64, b: f64, decimal_places: u8) -> bool {
    let factor = 10.0f64.powi(decimal_places as i32);
    let a = (a * factor).round();
    let b = (b * factor).round();
    a == b
}

#[test]
fn test_soil_tem_on_root_growth() {
    assert_eq!(SoilTemOnRootGrowth(30.), 1.);
    assert!(approx_equal(SoilTemOnRootGrowth(28.), 0.9712, 5));
    assert!(approx_equal(SoilTemOnRootGrowth(26.), 0.9168, 5));
    assert!(approx_equal(SoilTemOnRootGrowth(24.), 0.8368, 5));
    assert!(approx_equal(SoilTemOnRootGrowth(22.), 0.7312, 5));
    assert!(approx_equal(SoilTemOnRootGrowth(20.), 0.6000, 5));
    assert!(approx_equal(SoilTemOnRootGrowth(18.), 0.4432, 5));
    assert!(approx_equal(SoilTemOnRootGrowth(16.), 0.2608, 5));
    assert!(approx_equal(SoilTemOnRootGrowth(14.), 0.0528, 5));
    assert_eq!(SoilTemOnRootGrowth(13.5), 0.);
}

#[test]
fn test_doy_to_date() {
    let s = unsafe { str::from_utf8(slice::from_raw_parts(DoyToDate(1, 2021), 11)).unwrap() };
    assert_eq!(s, "01-JAN-2021");
    let s = unsafe { str::from_utf8(slice::from_raw_parts(DoyToDate(259, 1984), 11)).unwrap() };
    assert_eq!(s, "15-SEP-1984");
    let s = unsafe { str::from_utf8(slice::from_raw_parts(DoyToDate(0, 2021), 11)).unwrap() };
    assert_eq!(s, " ".repeat(11));
}

#[test]
fn test_date_to_doy() {
    assert_eq!(DateToDoy("01-JAN-2021".as_ptr() as *const i8, 2021), 1);
    assert_eq!(DateToDoy("01-OCT-2020".as_ptr() as *const i8, 2021), 275);
}

#[test]
fn test_leap_year() {
    assert_eq!(LeapYear(1900), 0);
    assert_eq!(LeapYear(2000), 1);
    assert_eq!(LeapYear(2019), 0);
    assert_eq!(LeapYear(2020), 1);
    assert_eq!(LeapYear(2021), 0);
}
