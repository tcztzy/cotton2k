#![feature(panic_always_abort)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
use io::toml::read_profile;
use pbr::ProgressBar;
mod bindings;
mod io;
mod profile;
use crate::bindings::{bEnd, C2KApp, DayFinish, DayStart, Daynum, WriteStateVariables};
use crate::io::to_csv::*;

pub fn run(profile_path: &std::path::Path) {
    let mut app: C2KApp = unsafe { C2KApp::new() };
    read_profile(profile_path).expect("Error in read_profile");
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
