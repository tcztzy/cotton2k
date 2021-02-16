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
extern "C" fn SlabLoc(isd: i32, row_space: f64) -> u32
//     This function computes the layer (lsdr) or column (ksdr) where the emitter 
//  of drip irrigation, or the fertilizer side - dressing is located. It is called
//  from ReadAgriculturalInput().
//     The following input arguments are used:
//               isd = horizontal or vertical distance
//     The following global variables are referenced here:
//                dl, nk, nl, wk.
//
{
    // horizontal
    if row_space > 0. {
        // Define the column of this location
        for w in 0..20 {
            if self::soil::width(w, row_space) >= isd as f64 {
                return w;
            }
        }
    } else {
        // Define the layer of this location
        for l in 0..40 {
            if self::soil::depth(l) >= isd as f64 {
                return l;
            }
        }
    }
    0
}
