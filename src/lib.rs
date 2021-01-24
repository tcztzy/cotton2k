#![feature(total_cmp)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

#[cfg(test)]
mod tests;

mod climate;
mod fruit;
mod leaf;
mod physiology;
mod root;
mod soil;
mod stem;
mod util;


include!(concat!(env!("OUT_DIR"), "/bindings.rs"));