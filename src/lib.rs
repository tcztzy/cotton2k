#[no_mangle]
pub extern fn LeapYear(year: u64) -> u64 {
    let is_divisible = |n| { year % n == 0 };

    (is_divisible(4) && (!is_divisible(100) || is_divisible(400))) as u64
}
