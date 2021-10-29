#![feature(panic_always_abort)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
use std::ffi::CString;
use std::path::Path;
mod bindings;
mod de;
mod io;
mod profile;
use crate::bindings::C2KApp;
use crate::io::toml::read_profile;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        panic!("profile file path should be provided!");
    }
    let mut app: C2KApp = unsafe { C2KApp::new() };
    read_profile(Path::new(&args[1])).expect("Error in read_profile");
    let filename = CString::new(args[2].to_string()).expect("Error");
    unsafe {
        app.RunTheModel(filename.as_ptr());
    }
}
