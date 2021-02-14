#![feature(total_cmp)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

#[cfg(test)]
mod tests;

mod climate;
mod fruit;
mod io;
mod leaf;
mod phenology;
mod physiology;
mod abscission;
mod root;
mod soil;
mod stem;
mod util;

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

#[no_mangle]
extern "C" fn SlabLoc(isd: i32, index: i32, wk: *const f64, dl: *const f64) -> i32
//     This function computes the layer (lsdr) or column (ksdr) where the emitter 
//  of drip irrigation, or the fertilizer side - dressing is located. It is called
//  from ReadAgriculturalInput().
//     The following input arguments are used:
//               isd = horizontal or vertical distance
//               index = 1 if horizontal or = 2 if vertical 
//     The following global variables are referenced here:
//                dl, nk, nl, wk.
//
{
    let wk = unsafe { std::slice::from_raw_parts(wk, 20) };
    let dl = unsafe { std::slice::from_raw_parts(dl, 40) };

    // horizontal
    if index == 1 {
        // Define the column of this location
        let mut sumwk = 0.; // sum of soil column widths.
        for w in wk.iter().enumerate() {
            sumwk += w.1;
            if sumwk >= isd as f64 {
                return w.0 as i32;
            }
        }
    } else if index == 2 {
        // Define the layer of this location
        let mut sumdl = 0.; // sum of soil layer depths.
        for d in dl.iter().enumerate() {
            sumdl += d.1;
            if sumdl >= isd as f64 {
                return d.0 as i32;
            }
        }
    }
    0
}
