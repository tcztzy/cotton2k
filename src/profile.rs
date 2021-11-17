use crate::de::{from_isoformat, from_isoformat_option};
use chrono::NaiveDate;
use serde::Deserialize;
use std::path::PathBuf;

#[inline]
fn zero() -> f64 {
    0.
}

#[inline]
fn zero_i32() -> i32 {
    0
}

#[derive(Deserialize, Debug)]
pub struct Profile {
    pub name: Option<String>,
    #[serde(default = "zero_i32")]
    pub light_intercept_method: i32,
    pub latitude: f64,
    pub longitude: f64,
    pub elevation: f64,
    #[serde(deserialize_with = "from_isoformat")]
    pub start_date: NaiveDate,
    #[serde(deserialize_with = "from_isoformat")]
    pub stop_date: NaiveDate,
    #[serde(default)]
    #[serde(deserialize_with = "from_isoformat_option")]
    pub emerge_date: Option<NaiveDate>,
    #[serde(default)]
    #[serde(deserialize_with = "from_isoformat_option")]
    pub plant_date: Option<NaiveDate>,
    pub co2_enrichment: Option<CO2Enrichment>,
    pub mulch: Option<Mulch>,
    pub weather_path: PathBuf,
    pub site: Site,
    pub cultivar_parameters: Vec<f64>,
    pub row_space: f64,
    #[serde(default = "zero")]
    pub skip_row_width: f64,
    pub plants_per_meter: f64,
    pub agronomy_operations: Vec<AgronomyOperation>,
    pub light_intercept_parameters: Option<[f64; 20]>,
    pub soil_layers: [SoilLayer; 14],
    pub soil_hydraulic: SoilHydraulic,
    pub plant_maps: Option<Vec<PlantMap>>,
}

#[derive(Deserialize, Debug)]
pub struct CO2Enrichment {
    pub factor: f64,
    #[serde(deserialize_with = "from_isoformat")]
    pub start_date: NaiveDate,
    #[serde(deserialize_with = "from_isoformat")]
    pub stop_date: NaiveDate,
}

#[derive(Deserialize, Debug)]
pub enum MulchType {
    NoMulch,
    All,                  // plastic layer on all soil surface
    OneColumnLeftAtSide, // plastic layer on all soil surface except one column at each side of the plant row.
    TwoColumnsLeftAtSide, // plastic layer on all soil surface except two columns at each side of the plant row.
}

#[derive(Deserialize, Debug)]
pub struct Mulch {
    pub indicator: MulchType,
    pub sw_trans: f64,
    pub lw_trans: f64,
    #[serde(deserialize_with = "from_isoformat")]
    pub start_date: NaiveDate,
    #[serde(default)]
    #[serde(deserialize_with = "from_isoformat_option")]
    pub stop_date: Option<NaiveDate>,
}

#[derive(Deserialize, Debug)]
pub struct Site {
    pub average_wind_speed: Option<f64>,
    pub estimate_dew_point: (f64, f64),
    pub wind_blow_after_sunrise: f64,
    pub wind_max_after_noon: f64,
    pub wind_stop_after_sunset: f64,
    pub night_time_wind_factor: f64,
    pub cloud_type_correction_factor: f64,
    pub max_temperature_after_noon: f64,
    pub deep_soil_temperature: (f64, f64, f64),
    pub dew_point_range: (f64, f64, f64),
    pub albedo_range: (f64, f64),
}

#[derive(Deserialize, Debug)]
pub struct WeatherRecord {
    #[serde(deserialize_with = "from_isoformat")]
    pub date: NaiveDate,
    pub irradiation: f64,
    pub tmax: f64,
    pub tmin: f64,
    pub rain: f64,
    pub wind: Option<f64>,
    pub tdew: Option<f64>,
}

#[inline]
fn default_predict() -> bool {
    false
}

#[inline]
fn default_irrgation_method() -> IrrigationMethod {
    IrrigationMethod::Sprinkler
}

#[inline]
fn default_fertilization_method() -> FertilizationMethod {
    FertilizationMethod::Broadcast
}

#[derive(Deserialize, Debug)]
pub enum IrrigationMethod {
    Sprinkler = 0,
    Furrow = 1,
    Drip = 2,
}

#[derive(Deserialize, Debug)]
pub enum FertilizationMethod {
    Broadcast = 0,
    Sidedress = 1,
    Foliar = 2,
    Drip = 3,
}

#[derive(Deserialize, Debug)]
pub enum PixMethod {
    Banded = 0,
    Sprinkler = 1,
    Broadcast = 2,
}

#[derive(Deserialize, Debug)]
#[serde(tag = "type")]
pub enum AgronomyOperation {
    irrigation {
        #[serde(deserialize_with = "from_isoformat")]
        date: NaiveDate,
        amount: f64,
        #[serde(default = "default_predict")]
        predict: bool,
        #[serde(default = "default_irrgation_method")]
        method: IrrigationMethod,
        #[serde(default = "zero")]
        drip_x: f64,
        #[serde(default = "zero")]
        drip_y: f64,
        max_amount: Option<f64>,
        #[serde(default)]
        #[serde(deserialize_with = "from_isoformat_option")]
        stop_predict_date: Option<NaiveDate>,
    },
    fertilization {
        #[serde(deserialize_with = "from_isoformat")]
        date: NaiveDate,
        #[serde(default = "zero")]
        urea: f64,
        #[serde(default = "zero")]
        nitrate: f64,
        #[serde(default = "zero")]
        ammonium: f64,
        #[serde(default = "default_fertilization_method")]
        method: FertilizationMethod,
        #[serde(default = "zero")]
        drip_x: f64,
        #[serde(default = "zero")]
        drip_y: f64,
    },
    defoliation {
        #[serde(deserialize_with = "from_isoformat")]
        date: NaiveDate,
        open_ratio: i32,
        #[serde(default = "default_predict")]
        predict: bool,
        #[serde(default = "zero")]
        ppa: f64, // pints per acre
    },
    cultivation {
        #[serde(deserialize_with = "from_isoformat")]
        date: NaiveDate,
        depth: f64,
    },
    pix {
        #[serde(deserialize_with = "from_isoformat")]
        date: NaiveDate,
        method: FertilizationMethod,
        ppa: f64, // pints per acre
    },
    watertable {
        #[serde(deserialize_with = "from_isoformat")]
        date: NaiveDate,
        level: f64,
        ecs: f64,
    },
}

#[derive(Deserialize, Debug)]
pub struct SoilLayer {
    pub ammonium: f64,
    pub nitrate: f64,
    pub organic_matter: f64,
    pub water_content: f64,
}

#[derive(Deserialize, Debug)]
pub struct SoilHydraulic {
    pub implicit_ratio: f64,
    pub max_conductivity: f64,
    pub psi_fc: f64,
    pub psi_id: f64,
    pub layers: Vec<SoilHydraulicLayer>,
}

#[derive(Deserialize, Debug)]
pub struct SoilHydraulicLayer {
    pub depth: f64,
    pub theta_d: f64,
    pub theta_s: f64,
    pub alpha: f64,
    pub beta: f64,
    pub hcs: f64,
    pub hcfc: f64,
    pub bulk_density: f64,
    pub clay: f64,
    pub sand: f64,
}

#[derive(Deserialize, Debug)]
pub struct PlantMap {
    #[serde(deserialize_with = "from_isoformat")]
    pub date: NaiveDate,
    pub plant_height: f64,
    pub main_stem_nodes: f64,
    pub number_of_squares: f64,
    pub number_of_bolls: f64,
    pub number_of_nodes: f64,
}
