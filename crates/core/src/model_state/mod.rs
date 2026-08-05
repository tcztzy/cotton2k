//! Explicit model-owned state grouped by simulation concern.
//!
//! [`ModelState`] contains calendar, plant, soil, atmosphere, flux, and
//! management state plus the legacy mirror used by translated routines. The
//! default state allocates the fixed Cotton2K grid sizes; the runner transfers
//! ownership of this value between daily steps and synchronizes the legacy
//! mirror only at domain boundaries.

use ndarray::{Array1, Array2, Array3};
mod grid_ops;
#[allow(missing_docs)]
mod legacy;
pub use grid_ops::{
    for_each_cell, for_each_cell_in_rows, for_each_fruiting_site, for_each_layer_col_span,
    for_each_row, for_each_row_mut,
};
pub use legacy::LegacyGlobalState;

/// Maximum active soil layers in the fixed simulation grid.
pub const MAX_LAYERS: usize = 40;
/// Maximum active soil columns in the fixed simulation grid.
pub const MAX_COLS: usize = 20;
/// Maximum vegetative branches represented by the legacy model.
pub const MAX_VEG_BRANCHES: usize = 3;
/// Maximum fruiting branches per vegetative branch.
pub const MAX_FRUIT_BRANCHES: usize = 30;
/// Maximum nodes per fruiting branch.
pub const MAX_NODES: usize = 5;
/// Maximum pre-fruiting nodes on the main stem.
pub const MAX_PREFRUIT_NODES: usize = 9;
/// Number of canopy light-interception layers.
pub const MAX_LAI_LAYERS: usize = 20;
/// Number of hourly atmospheric slots in one day.
pub const HOURS_PER_DAY: usize = 24;

#[derive(Debug, Clone)]
/// Calendar counters shared across daily simulation phases.
pub struct SimCalendarState {
    /// Current day of year.
    pub daynum: i32,
    /// First simulated day of year.
    pub day_start: i32,
    /// Last simulated day of year.
    pub day_finish: i32,
    /// Emergence day of year.
    pub day_emerge: i32,
    /// Days since emergence.
    pub kday: i32,
    /// Days since simulation start.
    pub day_of_simulation: i32,
}

impl Default for SimCalendarState {
    fn default() -> Self {
        Self {
            daynum: 0,
            day_start: 0,
            day_finish: 0,
            day_emerge: 0,
            kday: 0,
            day_of_simulation: 0,
        }
    }
}

#[derive(Debug, Clone)]
/// Plant-specific state extracted from the legacy plant arrays.
pub struct PlantState {
    /// Day of first square formation.
    pub first_square: i32,
    /// Number of pre-fruiting nodes currently active.
    pub num_pre_fru_nodes: i32,
    /// Number of vegetative branches currently active.
    pub num_veg_branches: i32,
    /// Fruiting-branch counts indexed by vegetative branch.
    pub num_fruit_branches: Array1<i32>,
    /// Node counts indexed by vegetative and fruiting branch.
    pub num_nodes: Array2<i32>,
    /// Pre-fruiting leaf weights.
    pub leaf_weight_pre_fru: Array1<f64>,
    /// Main-stem leaf weights by branch and site.
    pub leaf_weight_main_stem: Array2<f64>,
    /// Node leaf weights by branch, site, and node.
    pub leaf_weight_nodes: Array3<f64>,
    /// Pre-fruiting leaf areas.
    pub leaf_area_pre_fru: Array1<f64>,
    /// Main-stem leaf areas by branch and site.
    pub leaf_area_main_stem: Array2<f64>,
    /// Node leaf areas by branch, site, and node.
    pub leaf_area_nodes: Array3<f64>,
}

impl Default for PlantState {
    fn default() -> Self {
        Self {
            first_square: 0,
            num_pre_fru_nodes: 0,
            num_veg_branches: 0,
            num_fruit_branches: Array1::zeros(MAX_VEG_BRANCHES),
            num_nodes: Array2::zeros((MAX_VEG_BRANCHES, MAX_FRUIT_BRANCHES)),
            leaf_weight_pre_fru: Array1::zeros(MAX_PREFRUIT_NODES),
            leaf_weight_main_stem: Array2::zeros((MAX_VEG_BRANCHES, MAX_FRUIT_BRANCHES)),
            leaf_weight_nodes: Array3::zeros((MAX_VEG_BRANCHES, MAX_FRUIT_BRANCHES, MAX_NODES)),
            leaf_area_pre_fru: Array1::zeros(MAX_PREFRUIT_NODES),
            leaf_area_main_stem: Array2::zeros((MAX_VEG_BRANCHES, MAX_FRUIT_BRANCHES)),
            leaf_area_nodes: Array3::zeros((MAX_VEG_BRANCHES, MAX_FRUIT_BRANCHES, MAX_NODES)),
        }
    }
}

#[derive(Debug, Clone)]
/// Soil water, nitrogen, potential, and root-weight grid state.
pub struct SoilState {
    /// Volumetric water content by layer and column.
    pub vol_water_content: Array2<f64>,
    /// Volumetric nitrate content by layer and column.
    pub vol_no3n_content: Array2<f64>,
    /// Volumetric ammonium content by layer and column.
    pub vol_nh4n_content: Array2<f64>,
    /// Volumetric urea content by layer and column.
    pub vol_urean_content: Array2<f64>,
    /// Soil water potential by layer and column.
    pub soil_psi: Array2<f64>,
    /// Root weight by layer, column, and root age group.
    pub root_weight: Array3<f64>,
}

impl Default for SoilState {
    fn default() -> Self {
        Self {
            vol_water_content: Array2::zeros((MAX_LAYERS, MAX_COLS)),
            vol_no3n_content: Array2::zeros((MAX_LAYERS, MAX_COLS)),
            vol_nh4n_content: Array2::zeros((MAX_LAYERS, MAX_COLS)),
            vol_urean_content: Array2::zeros((MAX_LAYERS, MAX_COLS)),
            soil_psi: Array2::zeros((MAX_LAYERS, MAX_COLS)),
            root_weight: Array3::zeros((MAX_LAYERS, MAX_COLS, 3)),
        }
    }
}

#[derive(Debug, Clone)]
/// Hourly atmospheric forcing arrays for one simulation day.
pub struct AtmosphereState {
    /// Air temperature by hour.
    pub air_temp: Array1<f64>,
    /// Incoming radiation by hour.
    pub radiation: Array1<f64>,
    /// Relative humidity by hour.
    pub relative_humidity: Array1<f64>,
    /// Wind speed by hour.
    pub wind_speed: Array1<f64>,
}

impl Default for AtmosphereState {
    fn default() -> Self {
        Self {
            air_temp: Array1::zeros(HOURS_PER_DAY),
            radiation: Array1::zeros(HOURS_PER_DAY),
            relative_humidity: Array1::zeros(HOURS_PER_DAY),
            wind_speed: Array1::zeros(HOURS_PER_DAY),
        }
    }
}

#[derive(Debug, Clone)]
/// Daily fluxes and cumulative carbon/water exchange values.
pub struct FluxState {
    /// Net photosynthesis for the current day.
    pub net_photosynthesis: f64,
    /// Actual plant transpiration for the current day.
    pub actual_transpiration: f64,
    /// Actual soil evaporation for the current day.
    pub actual_soil_evaporation: f64,
    /// Cumulative net photosynthesis.
    pub cum_net_photosynth: f64,
}

impl Default for FluxState {
    fn default() -> Self {
        Self {
            net_photosynthesis: 0.0,
            actual_transpiration: 0.0,
            actual_soil_evaporation: 0.0,
            cum_net_photosynth: 0.0,
        }
    }
}

#[derive(Debug, Clone)]
/// Irrigation counters and predictive-irrigation settings.
pub struct ManagementState {
    /// Number of explicit irrigation operations.
    pub num_irrigations: i32,
    /// First day of predictive irrigation.
    pub day_start_pred_irrig: i32,
    /// Last day of predictive irrigation.
    pub day_stop_pred_irrig: i32,
    /// Encoded irrigation method.
    pub irrig_method: i32,
}

impl Default for ManagementState {
    fn default() -> Self {
        Self {
            num_irrigations: 0,
            day_start_pred_irrig: 0,
            day_stop_pred_irrig: 0,
            irrig_method: 0,
        }
    }
}

#[derive(Debug, Clone, Default)]
/// Complete explicit state passed between daily domain phases.
pub struct ModelState {
    /// Calendar counters.
    pub calendar: SimCalendarState,
    /// Plant arrays represented in explicit form.
    pub plant: PlantState,
    /// Soil arrays represented in explicit form.
    pub soil: SoilState,
    /// Atmospheric arrays represented in explicit form.
    pub atmosphere: AtmosphereState,
    /// Flux outputs and cumulative exchange.
    pub flux: FluxState,
    /// Management counters and settings.
    pub management: ManagementState,
    /// Compatibility mirror for translated legacy routines.
    pub legacy: LegacyGlobalState,
}

impl ModelState {
    /// Creates zeroed fixed-size state for a new simulation.
    pub fn new() -> Self {
        Self::default()
    }

    /// Sums leaf mass across pre-fruiting and fruiting sites.
    pub fn total_leaf_weight(&self) -> f64 {
        let mut result = 0.0;
        if self.plant.first_square <= 0 {
            result += 0.2;
        }

        for i in 0..self.plant.num_pre_fru_nodes.max(0) as usize {
            result += self.plant.leaf_weight_pre_fru[i];
        }

        for_each_fruiting_site(
            self.plant.num_veg_branches,
            &self.plant.num_fruit_branches,
            &self.plant.num_nodes,
            |k, l, nodes| {
                result += self.plant.leaf_weight_main_stem[[k, l]];
                result += self
                    .plant
                    .leaf_weight_nodes
                    .slice(ndarray::s![k, l, 0..nodes])
                    .sum();
            },
        );
        result
    }

    /// Sums leaf area across pre-fruiting and fruiting sites.
    pub fn total_leaf_area(&self) -> f64 {
        let mut result = 0.0;
        if self.plant.first_square <= 0 {
            result += 0.20 * 0.6;
        }

        for i in 0..self.plant.num_pre_fru_nodes.max(0) as usize {
            result += self.plant.leaf_area_pre_fru[i];
        }

        for_each_fruiting_site(
            self.plant.num_veg_branches,
            &self.plant.num_fruit_branches,
            &self.plant.num_nodes,
            |k, l, nodes| {
                result += self.plant.leaf_area_main_stem[[k, l]];
                result += self
                    .plant
                    .leaf_area_nodes
                    .slice(ndarray::s![k, l, 0..nodes])
                    .sum();
            },
        );
        result
    }
}
