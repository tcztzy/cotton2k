use ndarray::{Array1, Array2, Array3};
mod grid_ops;
mod legacy;
pub(crate) use grid_ops::{
    for_each_cell, for_each_cell_in_rows, for_each_fruiting_site, for_each_layer_col_span,
    for_each_row, for_each_row_mut,
};
pub use legacy::LegacyGlobalState;

pub const MAX_LAYERS: usize = 40;
pub const MAX_COLS: usize = 20;
pub const MAX_VEG_BRANCHES: usize = 3;
pub const MAX_FRUIT_BRANCHES: usize = 30;
pub const MAX_NODES: usize = 5;
pub const MAX_PREFRUIT_NODES: usize = 9;
pub const MAX_LAI_LAYERS: usize = 20;
pub const HOURS_PER_DAY: usize = 24;

#[derive(Debug, Clone)]
pub struct SimCalendarState {
    pub daynum: i32,
    pub day_start: i32,
    pub day_finish: i32,
    pub day_emerge: i32,
    pub kday: i32,
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
pub struct PlantState {
    pub first_square: i32,
    pub num_pre_fru_nodes: i32,
    pub num_veg_branches: i32,
    pub num_fruit_branches: Array1<i32>,
    pub num_nodes: Array2<i32>,
    pub leaf_weight_pre_fru: Array1<f64>,
    pub leaf_weight_main_stem: Array2<f64>,
    pub leaf_weight_nodes: Array3<f64>,
    pub leaf_area_pre_fru: Array1<f64>,
    pub leaf_area_main_stem: Array2<f64>,
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
pub struct SoilState {
    pub vol_water_content: Array2<f64>,
    pub vol_no3n_content: Array2<f64>,
    pub vol_nh4n_content: Array2<f64>,
    pub vol_urean_content: Array2<f64>,
    pub soil_psi: Array2<f64>,
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
pub struct AtmosphereState {
    pub air_temp: Array1<f64>,
    pub radiation: Array1<f64>,
    pub relative_humidity: Array1<f64>,
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
pub struct FluxState {
    pub net_photosynthesis: f64,
    pub actual_transpiration: f64,
    pub actual_soil_evaporation: f64,
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
pub struct ManagementState {
    pub num_irrigations: i32,
    pub day_start_pred_irrig: i32,
    pub day_stop_pred_irrig: i32,
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
pub struct ModelState {
    pub calendar: SimCalendarState,
    pub plant: PlantState,
    pub soil: SoilState,
    pub atmosphere: AtmosphereState,
    pub flux: FluxState,
    pub management: ManagementState,
    pub legacy: LegacyGlobalState,
}

impl ModelState {
    pub fn new() -> Self {
        Self::default()
    }

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
