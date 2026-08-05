//! C-compatible constants and records shared by the translated model boundary.
//!
//! The climate metric identifiers and fixed-layout records in this module are
//! consumed by the core global-state layer. They define representation and
//! array-size compatibility; they do not perform FFI calls or own storage.

use std::os::raw::{c_int, c_uint};

/// Maximum number of soil layers in the fixed legacy grid.
pub const maxl: c_int = 40;
/// Maximum number of soil columns in the fixed legacy grid.
pub const maxk: c_int = 20;
/// Legacy mathematical constant retained for translated equations.
pub const pi: f64 = 3.14159;

/// Numeric identifier type for climate-buffer fields.
pub type CLIMATE_METRIC = c_uint;
/// Daily maximum temperature metric.
pub const CLIMATE_METRIC_TMAX: CLIMATE_METRIC = 0;
/// Daily minimum temperature metric.
pub const CLIMATE_METRIC_TMIN: CLIMATE_METRIC = 1;
/// Daily irradiation metric.
pub const CLIMATE_METRIC_IRRD: CLIMATE_METRIC = 2;
/// Daily rainfall metric.
pub const CLIMATE_METRIC_RAIN: CLIMATE_METRIC = 3;
/// Daily wind metric.
pub const CLIMATE_METRIC_WIND: CLIMATE_METRIC = 4;
/// Daily dew-point temperature metric.
pub const CLIMATE_METRIC_TDEW: CLIMATE_METRIC = 5;

#[repr(C)]
#[derive(Debug, Copy, Clone, Default)]
/// Fixed-layout daily climate record used by the legacy buffer.
pub struct Climstruct {
    /// Day-of-year key.
    pub nDay: c_int,
    /// Irradiation in the legacy radiation unit.
    pub Rad: f64,
    /// Daily maximum air temperature.
    pub Tmax: f64,
    /// Daily minimum air temperature.
    pub Tmin: f64,
    /// Daily rainfall.
    pub Rain: f64,
    /// Daily wind value.
    pub Wind: f64,
    /// Daily dew-point temperature.
    pub Tdew: f64,
}

#[repr(C)]
#[derive(Debug, Copy, Clone, Default)]
/// Fixed-layout irrigation operation used by the legacy buffer.
pub struct Irrigation {
    /// Day-of-year of application.
    pub day: c_int,
    /// Encoded irrigation method.
    pub method: c_int,
    /// Drip column location.
    pub LocationColumnDrip: c_int,
    /// Drip layer location.
    pub LocationLayerDrip: c_int,
    /// Applied water amount in model units.
    pub amount: f64,
}
