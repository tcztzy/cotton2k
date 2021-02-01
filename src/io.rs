use std::ffi::CStr;
use std::fs::File;
use std::io::prelude::*;
use std::os::raw::c_char;
use chrono::prelude::*;


#[no_mangle]
extern "C" fn b01(profile_name: *const c_char, description: *const c_char) {
    let local: DateTime<Local> = Local::now();
    let profile_name = unsafe { CStr::from_ptr(profile_name).to_str().unwrap() };
    let description = unsafe { CStr::from_ptr(description).to_str().unwrap() };
    let mut file = File::create(format!("output/{}.B01", profile_name)).unwrap();
    let mut w = Vec::new();
    writeln!(&mut w, "{:>50}", "COTTON2K Version 4.0 (2003)").unwrap();
    writeln!(&mut w, "{:>62}", "A simulation model for irrigated cotton in arid regions").unwrap();
    writeln!(&mut w, "{:>50}", "Written by Avishalom Marani").unwrap();
    writeln!(&mut w).unwrap();
    writeln!(&mut w, "Profile Name:    {:<20}", profile_name).unwrap();
    writeln!(&mut w, "Simulation Date: {:<30}", local.format("%A, %B %d, %Y").to_string()).unwrap();
    writeln!(&mut w, "Description:     {:<55}", description).unwrap();
    writeln!(&mut w).unwrap();
    file.write_all(&mut w).unwrap();
}
