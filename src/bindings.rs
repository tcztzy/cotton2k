use std::os::raw::{c_int, c_uint};

pub const maxl: c_int = 40;
pub const maxk: c_int = 20;
pub const pi: f64 = 3.14159;

pub type CLIMATE_METRIC = c_uint;
pub const CLIMATE_METRIC_TMAX: CLIMATE_METRIC = 0;
pub const CLIMATE_METRIC_TMIN: CLIMATE_METRIC = 1;
pub const CLIMATE_METRIC_IRRD: CLIMATE_METRIC = 2;
pub const CLIMATE_METRIC_RAIN: CLIMATE_METRIC = 3;
pub const CLIMATE_METRIC_WIND: CLIMATE_METRIC = 4;
pub const CLIMATE_METRIC_TDEW: CLIMATE_METRIC = 5;

#[repr(C)]
#[derive(Debug, Copy, Clone, Default)]
pub struct Climstruct {
    pub nDay: c_int,
    pub Rad: f64,
    pub Tmax: f64,
    pub Tmin: f64,
    pub Rain: f64,
    pub Wind: f64,
    pub Tdew: f64,
}

#[repr(C)]
#[derive(Debug, Copy, Clone, Default)]
pub struct Irrigation {
    pub day: c_int,
    pub method: c_int,
    pub LocationColumnDrip: c_int,
    pub LocationLayerDrip: c_int,
    pub amount: f64,
}
