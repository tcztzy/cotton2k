#![feature(panic_always_abort)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
use pbr::ProgressBar;
use std::path::Path;
mod bindings;
mod io;
mod profile;
use crate::bindings::{bEnd, C2KApp, DayFinish, DayStart, Daynum, WriteStateVariables};
use crate::io::to_csv::*;
use crate::io::toml::read_profile;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        panic!("profile file path should be provided!");
    }
    let mut app: C2KApp = unsafe { C2KApp::new() };
    read_profile(Path::new(&args[1])).expect("Error in read_profile");
    output_file_headers().unwrap();
    unsafe {
        let count = (DayFinish - DayStart + 1) as u64;
        let mut pb = ProgressBar::new(count);
        // Do daily simulations
        Daynum = DayStart - 1;
        bEnd = false;
        //     Start the daily loop. If variable bEnd has been assigned a value
        //  of true end simulation.
        while !bEnd {
            let bAdjustToDo = app.DoAdjustments();
            pb.inc();
            //     Execute simulation for this day.
            app.SimulateThisDay();
            write_record().unwrap();
            //     If there are pending plant adjustments, call
            //     WriteStateVariables() to write
            //  state variables of this day in a scratch file.
            if bAdjustToDo {
                WriteStateVariables(true);
            }
        } // end while
    }
}
