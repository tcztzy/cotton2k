#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
use std::ffi::CString;
use std::path::Path;
use std::io::Read;
use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Debug)]
struct Profile {
    latitude: f64,
    longitude: f64,
    elevation: f64,
}

fn read_profile(profile_path: &Path) {
    let mut file = std::fs::File::open(profile_path).unwrap();
    let mut contents = String::new();
    file.read_to_string(&mut contents).unwrap();
    let profile: Profile = toml::from_str(&contents).unwrap();
    unsafe {
        Latitude = profile.latitude;
        Longitude = profile.longitude;
        Elevation = profile.elevation;
    }
}


include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        panic!("profile file path should be provided!");
    }
    let mut app: C2KApp = unsafe { C2KApp::new() };
    read_profile(Path::new(&args[1]));
    let filename = CString::new(args[2].to_string()).expect("Error");
    unsafe {
        app.RunTheModel(filename.as_ptr());
    }
}
