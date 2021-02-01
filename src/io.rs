use chrono::prelude::*;
use std::ffi::CStr;
use std::fs::{File, OpenOptions};
use std::io::prelude::*;
use std::os::raw::c_char;

#[no_mangle]
extern "C" fn b01(profile_name: *const c_char, description: *const c_char) {
    let local: DateTime<Local> = Local::now();
    let profile_name = unsafe { CStr::from_ptr(profile_name).to_str().unwrap() };
    let description = unsafe { CStr::from_ptr(description).to_str().unwrap() };
    let mut file = File::create(format!("output/{}.B01", profile_name)).unwrap();
    writeln!(file, "{:>50}", "COTTON2K Version 4.0 (2003)").unwrap();
    writeln!(
        file,
        "{:>62}",
        "A simulation model for irrigated cotton in arid regions"
    )
    .unwrap();
    writeln!(file, "{:>50}", "Written by Avishalom Marani").unwrap();
    writeln!(file).unwrap();
    writeln!(file, "Profile Name:    {:<20}", profile_name).unwrap();
    writeln!(
        file,
        "Simulation Date: {:<30}",
        local.format("%A, %B %d, %Y").to_string()
    )
    .unwrap();
    writeln!(file, "Description:     {:<55}", description).unwrap();
    writeln!(file).unwrap();
}

#[no_mangle]
extern "C" fn b01_append_defoliant(profile_name: *const c_char, date: *const c_char) {
    let profile_name = unsafe { CStr::from_ptr(profile_name).to_str().unwrap() };
    let date = unsafe { CStr::from_ptr(date).to_str().unwrap() };
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(format!("output/{}.B01", profile_name))
        .unwrap();
    writeln!(file, " ****   defoliant applied on {}    ****", date).unwrap();
}
