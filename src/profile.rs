use crate::input_functions::form;
use crate::plant::root::RootImpedanceTables;
use crate::*;
use chrono::Datelike;
use serde::Deserialize;
use std::io::Write;
use std::path::PathBuf;

struct NaiveDateVisitor;

impl<'de> serde::de::Visitor<'de> for NaiveDateVisitor {
    type Value = chrono::NaiveDate;

    fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(formatter, "a string represents chrono::NaiveDate")
    }

    fn visit_str<E>(self, s: &str) -> Result<Self::Value, E>
    where
        E: serde::de::Error,
    {
        match chrono::NaiveDate::parse_from_str(s, "%F") {
            Ok(t) => Ok(t),
            Err(_) => Err(serde::de::Error::invalid_value(
                serde::de::Unexpected::Str(s),
                &self,
            )),
        }
    }
}
struct OptionalNaiveDateVisitor;

impl<'de> serde::de::Visitor<'de> for OptionalNaiveDateVisitor {
    type Value = Option<chrono::NaiveDate>;

    fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(formatter, "a optional string represents chrono::NaiveDate")
    }

    fn visit_none<E>(self) -> Result<Self::Value, E>
    where
        E: serde::de::Error,
    {
        Ok(None)
    }

    fn visit_some<D>(self, deserializer: D) -> Result<Self::Value, D::Error>
    where
        D: serde::de::Deserializer<'de>,
    {
        Ok(Some(deserializer.deserialize_str(NaiveDateVisitor)?))
    }
}

fn from_isoformat<'de, D>(d: D) -> Result<chrono::NaiveDate, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    d.deserialize_str(NaiveDateVisitor)
}

fn from_isoformat_option<'de, D>(d: D) -> Result<Option<chrono::NaiveDate>, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    d.deserialize_option(OptionalNaiveDateVisitor)
}

#[inline]
fn zero() -> f64 {
    0.
}

#[inline]
fn original_light_intercept_method() -> LightInterceptMethod {
    LightInterceptMethod::Original
}

#[inline]
fn default_last_day_weather_data() -> NaiveDate {
    NaiveDate::from_yo_opt(1970, 1).unwrap()
}

#[derive(Deserialize, Debug)]
pub struct Profile {
    #[serde(skip)]
    pub path: PathBuf,
    #[serde(default = "original_light_intercept_method")]
    pub light_intercept_method: LightInterceptMethod,
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
    pub soil_impedance: Option<PathBuf>,
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
    /// maximum leaf area index.
    #[serde(skip)]
    pub lmax: f64,
    /// The effect of moisture stress on the photosynthetic rate
    #[serde(skip)]
    pub ptsred: f64,
    /// correction factor for ambient CO2 in air
    #[serde(skip)]
    pub ambient_CO2_factor: f64,
    #[serde(skip)]
    pub num_watertable_data: usize,
    #[serde(skip)]
    pub states: Vec<State>,
    #[serde(skip)]
    #[serde(default = "default_last_day_weather_data")]
    pub last_day_weather_data: NaiveDate,
    #[serde(skip)]
    pub root_impedance_tables: Option<RootImpedanceTables>,
    #[serde(skip)]
    pub model_state: ModelState,
}

#[derive(Deserialize, Debug, Clone, Copy)]
pub struct CO2Enrichment {
    pub factor: f64,
    #[serde(deserialize_with = "from_isoformat")]
    pub start_date: NaiveDate,
    #[serde(deserialize_with = "from_isoformat")]
    pub stop_date: NaiveDate,
}

impl Default for CO2Enrichment {
    fn default() -> Self {
        CO2Enrichment {
            factor: 1.,
            start_date: NaiveDate::from_ymd_opt(1900, 1, 1).unwrap(),
            stop_date: NaiveDate::from_ymd_opt(2100, 1, 1).unwrap(),
        }
    }
}

#[derive(Deserialize, Debug, Clone, Copy)]
pub enum MulchType {
    NoMulch,
    /// plastic layer on all soil surface
    All,
    /// plastic layer on all soil surface except one column at each side of the plant row.
    OneColumnLeftAtSide,
    /// plastic layer on all soil surface except two columns at each side of the plant row.
    TwoColumnsLeftAtSide,
}

#[derive(Deserialize, Debug, Clone, Copy)]
pub enum LightInterceptMethod {
    Original,
    Fry1980,
    Latered,
}

#[derive(Deserialize, Debug, Clone, Copy)]
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
fn default_irrigation_method() -> IrrigationMethod {
    IrrigationMethod::Sprinkler
}

#[inline]
fn default_fertilization_method() -> FertilizationMethod {
    FertilizationMethod::Broadcast
}

#[derive(Deserialize, Debug, Clone, Copy)]
pub enum IrrigationMethod {
    Sprinkler = 0,
    Furrow = 1,
    Drip = 2,
}

#[derive(Deserialize, Debug, Clone, Copy)]
pub enum FertilizationMethod {
    Broadcast = 0,
    Sidedress = 1,
    Foliar = 2,
    Drip = 3,
}

#[derive(Deserialize, Debug, Clone, Copy)]
#[serde(tag = "type")]
pub enum AgronomyOperation {
    irrigation(AgronomyOperationIrrigation),
    /// nitrogen fertilizer application information for each day.
    fertilization {
        /// date of application
        #[serde(deserialize_with = "from_isoformat")]
        date: NaiveDate,
        /// urea N applied, kg N per ha
        #[serde(default = "zero")]
        urea: f64,
        /// nitrate N applied, kg N per ha;
        #[serde(default = "zero")]
        nitrate: f64,
        /// ammonium N applied, kg N per ha;
        #[serde(default = "zero")]
        ammonium: f64,
        /// method of application
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
        /// pints per acre
        #[serde(default = "zero")]
        ppa: f64,
    },
    cultivation {
        #[serde(deserialize_with = "from_isoformat")]
        date: NaiveDate,
        depth: f64,
    },
    watertable(AgronomyOperationWaterTable),
}
#[derive(Deserialize, Debug, Clone, Copy)]
pub struct AgronomyOperationIrrigation {
    #[serde(deserialize_with = "from_isoformat")]
    date: NaiveDate,
    amount: f64,
    #[serde(default = "default_predict")]
    predict: bool,
    #[serde(default = "default_irrigation_method")]
    method: IrrigationMethod,
    #[serde(default = "zero")]
    drip_x: f64,
    #[serde(default = "zero")]
    drip_y: f64,
    max_amount: Option<f64>,
    #[serde(default)]
    #[serde(deserialize_with = "from_isoformat_option")]
    stop_predict_date: Option<NaiveDate>,
}

#[derive(Deserialize, Debug, Clone, Copy)]
pub struct AgronomyOperationWaterTable {
    #[serde(deserialize_with = "from_isoformat")]
    pub date: NaiveDate,
    /// water table level input data (cm below soil surface).
    pub level: f64,
    /// electrical conductivity of saturated soil extract (mmho/cm)
    pub ecs: f64,
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

/// estimates the approximate daily average dewpoint temperature when it is not available.
fn estimate_dew_point(maxt: f64, site5: f64, site6: f64) -> f64 {
    if maxt <= 20. {
        site5
    } else if maxt >= 40. {
        site6
    } else {
        ((40. - maxt) * site5 + (maxt - 20.) * site6) / 20.
    }
}

/// This function opens the initial soil data file and reads it.
/// It is executed once at the beginning of the simulation.
/// It is called by [Profile::initialize()].
fn init_soil(
    legacy: &mut LegacyGlobalState,
    soil_layers: &[SoilLayer; 14],
    soil_hydraulic: &SoilHydraulic,
) {
    let mut condfc = [0f64; 9]; // hydraulic conductivity at field capacity of horizon layers, cm per day.
    let mut h2oint = [0f64; 14]; // initial soil water content, percent of field capacity,
                                 // defined by input for consecutive 15 cm soil layers.
    let mut ldepth = [0f64; 9]; // depth from soil surface to the end of horizon layers, cm.
    let mut oma = [0f64; 14]; // organic matter at the beginning of the season, percent of
                              // soil weight, defined by input for consecutive 15 cm soil
                              // layers.
    let mut pclay = [0f64; 9]; // percentage of clay in soil horizon of horizon layers.
    let mut psand = [0f64; 9]; // percentage of sand in soil horizon of horizon layers.
    let mut psidra: f64; // soil matric water potential, bars, for which immediate
                         // drainage will be simulated (suggested value -0.25 to -0.1).
    let mut psisfc: f64; // soil matric water potential at field capacity,
                         // bars (suggested value -0.33 to -0.1).
    let mut rnnh4 = [0f64; 14]; // residual nitrogen as ammonium in soil at beginning of
                                // season, kg per ha. defined by input for consecutive 15 cm
                                // soil layers.
    let mut rnno3 = [0f64; 14]; // residual nitrogen as nitrate in soil at beginning of
                                // season, kg per ha. defined by input for consecutive 15 cm
                                // soil layers.
    for (i, layer) in soil_layers.iter().enumerate() {
        rnnh4[i] = layer.ammonium;
        rnno3[i] = layer.nitrate;
        oma[i] = layer.organic_matter;
        h2oint[i] = layer.water_content;
    }
    let lyrsol = soil_hydraulic.layers.len();
    //     Zeroise arrays of data values.
    for i in 0..9 {
        legacy.airdr[i] = 0.;
        legacy.thetas[i] = 0.;
        legacy.alpha[i] = 0.;
        legacy.beta[i] = 0.;
        legacy.saturated_hyd_cond[i] = 0.;
        legacy.bulk_density[i] = 0.;
    }
    legacy.ratio_implicit = soil_hydraulic.implicit_ratio;
    legacy.conmax = soil_hydraulic.max_conductivity;
    psisfc = soil_hydraulic.psi_fc;
    psidra = soil_hydraulic.psi_id;
    if psisfc > 0. {
        psisfc = -psisfc; // make sure it is negative
    }
    if psidra > 0. {
        psidra = -psidra; // make sure it is negative
    }
    for (il, layer) in soil_hydraulic.layers.iter().enumerate() {
        ldepth[il] = layer.depth;
        condfc[il] = layer.hcfc;
        pclay[il] = layer.clay;
        psand[il] = layer.sand;
        legacy.airdr[il] = layer.theta_d;
        legacy.thetas[il] = layer.theta_s;
        legacy.alpha[il] = layer.alpha;
        legacy.beta[il] = layer.beta;
        legacy.saturated_hyd_cond[il] = layer.hcs;
        legacy.bulk_density[il] = layer.bulk_density;
    }

    let mut j = 0usize; // horizon number
    let mut sumdl = 0f64; // depth to the bottom this layer (cm);
    let rm = 2.65f64; // density of the solid fraction of the soil (g / cm3)
    let mut bdl = [0f64; 40]; // array of bulk density of soil layers

    for l in 0..legacy.nl as usize {
        //     Using the depth of each horizon layer (ldepth), the horizon
        //  number (legacy.soil_horizon_num) is computed for each soil layer.
        sumdl += legacy.dl[l];
        while (sumdl > ldepth[j]) && (j < lyrsol) {
            j += 1;
        }
        legacy.soil_horizon_num[l] = j as i32;
        // bdl, legacy.thad, legacy.thts are defined for each soil layer, using the respective input variables legacy.bulk_density, legacy.airdr,
        // legacy.thetas.
        //
        // legacy.field_capacity, legacy.max_water_capacity and legacy.thetar are computed for each layer, as water content ($cm^3\ cm^{-3}$)
        // of each layer corresponding to matric potentials of psisfc (for field capacity), psidra (for free drainage)
        // and -15 bars (for permanent wilting point), respectively, using function qpsi.
        //
        // pore space volume (legacy.pore_space) is also computed for each layer. make sure that saturated water content is not
        // more than pore space.
        bdl[l] = legacy.bulk_density[j];
        legacy.pore_space[l] = 1. - legacy.bulk_density[j] / rm;
        if legacy.thetas[j] > legacy.pore_space[l] {
            legacy.thetas[j] = legacy.pore_space[l];
        }
        legacy.thad[l] = legacy.airdr[j];
        legacy.thts[l] = legacy.thetas[j];
        legacy.field_capacity[l] = crate::general_functions::qpsi(
            psisfc,
            legacy.thad[l],
            legacy.thts[l],
            legacy.alpha[j],
            legacy.beta[j],
        );
        legacy.max_water_capacity[l] = crate::general_functions::qpsi(
            psidra,
            legacy.thad[l],
            legacy.thts[l],
            legacy.alpha[j],
            legacy.beta[j],
        );
        legacy.thetar[l] = crate::general_functions::qpsi(
            -15.,
            legacy.thad[l],
            legacy.thts[l],
            legacy.alpha[j],
            legacy.beta[j],
        );
        // When the saturated hydraulic conductivity (legacy.saturated_hyd_cond) is not given, it is computed from the hydraulic
        // conductivity at field capacity (condfc), using the wcond function.
        if legacy.saturated_hyd_cond[j] <= 0. {
            legacy.saturated_hyd_cond[j] = condfc[j]
                / crate::general_functions::wcond(
                    legacy.field_capacity[l],
                    legacy.thad[l],
                    legacy.thts[l],
                    legacy.beta[j],
                    1.,
                    1.,
                );
        }
    }
    // Loop for all soil layers. Compute depth from soil surface to the end of each layer (sumdl).
    sumdl = 0.;
    for l in 0..legacy.nl as usize {
        sumdl += legacy.dl[l];
        // At start of simulation compute estimated movable fraction of nitrates in each soil layer, following the work
        // of:
        //     Bowen, W.T., Jones, J.W., Carsky, R.J., and Quintana, J.O. 1993.
        //  Evaluation of the nitrogen submodel of CERES-maize following legume
        //  green manure incorporation. Agron. J. 85:153-159.
        //
        // The fraction of total nitrate in a layer that is in solution and can move from one layer to the next with
        // the downward flow of water is a function of the adsorption coefficient, soil bulk density, and the
        // volumetric soil water content at the drained upper limit.
        //
        // Adsorption coefficients are assumed to be 0.0 up to 30 cm depth, and deeper than 30 cm - 0.2, 0.4, 0.8, 1.0,
        // 1.2, and 1.6 for each successive 15 cm layer.

        let coeff: f64 = if sumdl <= 30. {
            0.
        } else if sumdl <= 45. {
            0.2
        } else if sumdl <= 60. {
            0.4
        } else if sumdl <= 75. {
            0.6
        } else if sumdl <= 90. {
            0.8
        } else if sumdl <= 105. {
            1.0
        } else if sumdl <= 120. {
            1.2
        } else {
            1.6
        };
        legacy.no3_flow_fraction[l] = 1. / (1. + coeff * bdl[l] / legacy.max_water_capacity[l]);
        // Determine the corresponding 15 cm layer of the input file. Compute the initial volumetric water content
        // (legacy.vol_water_content) of each layer, and check that it will not be less than the air-dry value or more than pore
        // space volume.
        j = ((sumdl - 1.) / 15.).floor() as usize;
        if j > 13 {
            j = 13;
        }
        let n = legacy.soil_horizon_num[l] as usize;
        legacy.vol_water_content[[l, 0]] = legacy.field_capacity[l] * h2oint[j] / 100.;
        if legacy.vol_water_content[[l, 0]] < legacy.airdr[n] {
            legacy.vol_water_content[[l, 0]] = legacy.airdr[n];
        }
        if legacy.vol_water_content[[l, 0]] > legacy.pore_space[l] {
            legacy.vol_water_content[[l, 0]] = legacy.pore_space[l];
        }
        // Initial values of ammonium N (rnnh4, legacy.vol_nh4_n_content) and nitrate N (rnno3, legacy.vol_no3_n_content) are converted
        // from kgs per ha to $mg\ cm^{-3}$ for each soil layer, after checking for minimal amounts.
        if rnno3[j] < 2.0 {
            rnno3[j] = 2.0;
        }
        if rnnh4[j] < 0.2 {
            rnnh4[j] = 0.2;
        }
        legacy.vol_no3_n_content[[l, 0]] = rnno3[j] / 15. * 0.01;
        legacy.vol_nh4_n_content[[l, 0]] = rnnh4[j] / 15. * 0.01;
        // organic matter in mg / cm3 units.
        let om = (oma[j] / 100.) * bdl[l] * 1000.;
        // potom is the proportion of readily mineralizable om. it is a function of soil depth (sumdl, in cm), modified
        // from GOSSYM (where it probably includes the 0.4 factor for organic C in om).
        let mut potom = 0.15125 - 0.02878 * sumdl.ln();
        if potom < 0. {
            potom = 0.;
        }
        // legacy.fresh_organic_matter is the readily mineralizable organic matter (="fresh organic matter" in CERES models).
        // legacy.humus_organic_matter is the remaining organic matter, which is mineralized very slowly.
        legacy.fresh_organic_matter[[l, 0]] = om * potom;
        legacy.humus_organic_matter[[l, 0]] = om * (1. - potom);
    }
    // Since the initial value has been set for the first column only in each layer, these values are now assigned to
    // all the other columns.
    for l in 0..legacy.nl as usize {
        for k in 1..legacy.nk as usize {
            legacy.vol_water_content[[l, k]] = legacy.vol_water_content[[l, 0]];
            legacy.vol_no3_n_content[[l, k]] = legacy.vol_no3_n_content[[l, 0]];
            legacy.vol_nh4_n_content[[l, k]] = legacy.vol_nh4_n_content[[l, 0]];
            legacy.fresh_organic_matter[[l, k]] = legacy.fresh_organic_matter[[l, 0]];
            legacy.humus_organic_matter[[l, k]] = legacy.humus_organic_matter[[l, 0]];
        }
    }
    // Total amounts of water (legacy.initial_total_soil_water), nitrate N (legacy.total_soil_no3_n), ammonium N (legacy.total_soil_nh4_n), and urea
    // N (legacy.total_soil_urea_n) are computed for the whole slab.
    legacy.initial_total_soil_water = 0.;
    legacy.total_soil_no3_n = 0.;
    legacy.total_soil_nh4_n = 0.;
    legacy.total_soil_urea_n = 0.;

    for l in 0..legacy.nl as usize {
        for k in 0..legacy.nk as usize {
            legacy.initial_total_soil_water +=
                legacy.vol_water_content[[l, k]] * legacy.dl[l] * legacy.wk[k];
            legacy.total_soil_no3_n +=
                legacy.vol_no3_n_content[[l, k]] * legacy.dl[l] * legacy.wk[k];
            legacy.total_soil_nh4_n +=
                legacy.vol_nh4_n_content[[l, k]] * legacy.dl[l] * legacy.wk[k];
            legacy.vol_urea_n_content[[l, k]] = 0.;
        }
    }
    // legacy.initial_total_soil_water is converted from cm3 per slab to mm.
    legacy.initial_total_soil_water = 10. * legacy.initial_total_soil_water / legacy.row_space;
    let bsand = 20f64; // heat conductivity of sand and silt (mcal cm-1 s-1 C-1).
    let bclay = 7f64; // heat conductivity of clay (mcal cm-1 s-1 C-1).
    let cka = 0.0615f64; // heat conductivity of air (mcal cm-1 s-1 C-1).
    let ckw = 1.45f64; // heat conductivity of water (mcal cm-1 s-1 C-1).
    let cmin = 0.46f64; // heat capacity of the mineral fraction of the soil.
    let corg = 0.6f64; // heat capacity of the organic fraction of the soil.
    let ga = 0.144f64; // shape factor for air in pore spaces.
    let ro = 1.3f64; // specific weight of organic fraction of soil.

    // Compute aggregation factors:
    legacy.dsand = form(bsand, ckw, ga); // aggregation factor for sand in water
    legacy.dclay = form(bclay, ckw, ga); // aggregation factor for clay in water
    let dsandair: f64 = form(bsand, cka, ga); // aggregation factor for sand in air
    let dclayair: f64 = form(bclay, cka, ga); // aggregation factor for clay in air

    // Loop over all soil layers, and define indices for some soil arrays.

    sumdl = 0.; // sum of depth of consecutive soil layers.

    for l in 0..legacy.nl as usize {
        sumdl += legacy.dl[l];
        let mut j = ((sumdl + 14.) / 15.).floor() as usize - 1; // layer definition for oma
        if j > 13 {
            j = 13;
        }
        // Using the values of the clay and organic matter percentages in the soil, compute mineral and organic
        // fractions of the soil, by weight and by volume.
        let mmo = oma[j] / 100.; // organic matter fraction of dry soil (by weight).
        let mm = 1. - mmo; // mineral fraction of dry soil (by weight).

        // legacy.marginal_water_content is set as a function of the sand fraction of the soil.
        let ra = (mmo / ro) / (mm / rm); // volume ratio of organic to mineral soil fractions.

        let i1 = legacy.soil_horizon_num[l] as usize; //  layer definition as in soil hydrology input file.

        // The volume fractions of clay (legacy.clay_volume_fraction) and of sand plus silt (legacy.sand_volume_fraction), are calculated
        legacy.marginal_water_content[l] = 0.1 - 0.07 * psand[i1] / 100.;
        let xo = (1. - legacy.pore_space[l]) * ra / (1. + ra); // organic fraction of soil (by volume).
        let xm = (1. - legacy.pore_space[l]) - xo; // mineral fraction of soil (by volume).
        legacy.clay_volume_fraction[l] = pclay[i1] * xm / mm / 100.;
        legacy.sand_volume_fraction[l] = 1. - legacy.pore_space[l] - legacy.clay_volume_fraction[l];
        // Heat capacity of the solid soil fractions (mineral + organic, by volume )
        legacy.heat_capacity_soil_solid[l] = xm * cmin + xo * corg;
        // The heat conductivity of dry soil (legacy.heat_cond_dry_soil) is computed using the procedure suggested by De Vries.
        legacy.heat_cond_dry_soil[l] = 1.25
            * (legacy.pore_space[l] * cka
                + dsandair * bsand * legacy.sand_volume_fraction[l]
                + dclayair * bclay * legacy.clay_volume_fraction[l])
            / (legacy.pore_space[l]
                + dsandair * legacy.sand_volume_fraction[l]
                + dclayair * legacy.clay_volume_fraction[l]);
    }
}

impl Profile {
    /// Run this profile.
    pub fn run(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        crate::runner::execute_profile(self, || false, |_, _| {})?;
        Ok(())
    }

    pub fn initialize(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        self.lmax = 0.;
        self.num_watertable_data = 0;
        self.model_state = ModelState::new();
        {
            let legacy = &mut self.model_state.legacy;
            initialize_global(legacy);
            legacy.iyear = self.start_date.year();
            legacy.day_start = self.start_date.ordinal() as i32;
            legacy.day_finish = self.stop_date.ordinal() as i32;
            match self.emerge_date {
                Some(date) => legacy.day_emerge = date.ordinal() as i32,
                None => {}
            }
            match self.plant_date {
                Some(date) => {
                    legacy.day_plant = date.ordinal() as i32;
                }
                None => {}
            }
            if self.emerge_date.is_none() {
                // If the date of emergence has not been given, emergence will be simulated by the model. In this case,
                // legacy.isw = 0, and a check is performed to make sure that the date of planting has been given.
                if self.plant_date.is_none() {
                    panic!(
                        "one of planting date or emergence date must be given in the profile file!!"
                    );
                }
                legacy.isw = 0;
            } else if self.emerge_date > self.plant_date {
                // If the date of emergence has been given in the input: legacy.isw = 1 if simulation starts before emergence,
                legacy.isw = 1;
            } else {
                // or legacy.isw = 2 if simulation starts at emergence.
                legacy.isw = 2;
                legacy.kday = 1;
            }
            // For advanced users only: if there is CO2 enrichment, read also CO2 factor, DOY dates
            // If soil mulch is used, read relevant parameters.
            match self.mulch {
                Some(mulch) => match mulch.indicator {
                    MulchType::NoMulch => {}
                    _ => {
                        legacy.mulch_indicator = mulch.indicator as i32;
                        legacy.mulch_tran_sw = mulch.sw_trans;
                        legacy.mulch_tran_lw = mulch.lw_trans;
                        legacy.day_start_mulch = mulch.start_date.ordinal() as i32;
                        match mulch.stop_date {
                            Some(date) => legacy.day_end_mulch = date.ordinal() as i32,
                            None => legacy.day_end_mulch = legacy.day_finish,
                        }
                    }
                },
                None => {
                    legacy.mulch_indicator = MulchType::NoMulch as i32;
                }
            }
            for pair in self.cultivar_parameters.iter().enumerate() {
                legacy.var_par[pair.0 + 1] = *pair.1;
            }
            legacy.row_space = if self.skip_row_width > 0. {
                (self.row_space + self.skip_row_width) / 2.
            } else {
                self.row_space
            };
            // legacy.plant_row_location is the distance from edge of slab, cm, of the plant row.
            legacy.plant_row_location = legacy.row_space / 2.;
            // Compute legacy.plant_population - number of plants per hectar, and legacy.per_plant_area - the average surface area per
            // plant, in $dm^2$, and the empirical plant density factor (legacy.density_factor). This factor will be used to
            // express the effect of plant density on some plant growth rate functions.
            // Note that legacy.density_factor =1 for 5 plants per sq m (or 50000 per ha).
            legacy.plant_population = self.plants_per_meter / legacy.row_space * 1000000.;
            legacy.per_plant_area = 1000000. / legacy.plant_population;
            legacy.density_factor =
                (legacy.var_par[1] * (5. - legacy.plant_population / 10000.)).exp();
            // Define the numbers of rows and columns in the soil slab (legacy.nl, legacy.nk).
            // Define the depth, in cm, of consecutive legacy.nl layers.
            legacy.nl = maxl;
            legacy.nk = maxk;
            legacy.dl[0] = 2.;
            legacy.dl[1] = 2.;
            legacy.dl[2] = 2.;
            legacy.dl[3] = 4.;
            for i in 4..(maxl - 2) as usize {
                legacy.dl[i] = 5.;
            }
            legacy.dl[(maxl - 2) as usize] = 10.;
            legacy.dl[(maxl - 1) as usize] = 10.;
            //      The width of the slab columns is computed by dividing the row
            //  spacing by the number of columns. It is assumed that slab width is
            //  equal to the average row spacing, and column widths are uniform.
            //      Note: legacy.wk is an array - to enable the option of non-uniform
            //  column widths in the future.
            //      legacy.plant_row_column (the column including the plant row) is now computed
            //      from
            //  legacy.plant_row_location (the distance of the plant row from the edge of the
            //  slab).
            let mut sumwk = 0.; // sum of column widths
            legacy.plant_row_column = 0;
            for k in 0..legacy.nk {
                legacy.wk[k as usize] = legacy.row_space / legacy.nk as f64;
                sumwk = sumwk + legacy.wk[k as usize];
                if legacy.plant_row_column == 0 && sumwk > legacy.plant_row_location {
                    legacy.plant_row_column =
                        if (sumwk - legacy.plant_row_location) > (0.5 * legacy.wk[k as usize]) {
                            k - 1
                        } else {
                            k
                        };
                }
            }
        }
        let mut rdr = csv::Reader::from_path(&self.weather_path)?;
        let mut jdd: u32 = 0;
        let day_start = self.start_date.ordinal() as i32;
        let sim_year = self.start_date.year();
        {
            let mut clim = Clim.write().expect("Clim lock poisoned");
            for result in rdr.deserialize() {
                let record: WeatherRecord = result?;
                jdd = record.date.ordinal();
                let j = jdd as i32 - day_start;
                if j < 0 {
                    continue;
                }
                clim[j as usize].nDay = jdd as i32;
                // convert $MJ\ m^{-2}$ to langleys
                clim[j as usize].Rad = record.irradiation * 23.884;
                clim[j as usize].Tmax = record.tmax;
                clim[j as usize].Tmin = record.tmin;
                clim[j as usize].Wind =
                    if self.site.average_wind_speed.is_some() && record.wind.is_none() {
                        self.site.average_wind_speed.unwrap()
                    } else {
                        record.wind.unwrap_or(0.)
                    };
                clim[j as usize].Tdew = record.tdew.unwrap_or(estimate_dew_point(
                    record.tmax,
                    self.site.estimate_dew_point.0,
                    self.site.estimate_dew_point.1,
                ));
                clim[j as usize].Rain = record.rain;
            }
        }
        self.last_day_weather_data = NaiveDate::from_yo_opt(sim_year, jdd).unwrap();
        let mut idef: usize = 0;
        let impedance_tables = self.read_soil_impedance(self.soil_impedance.as_ref().unwrap())?;
        self.root_impedance_tables = Some(impedance_tables);
        {
            let legacy = &mut self.model_state.legacy;
            legacy.num_irrigations = 0;
            for i in 0..5 {
                legacy.defoliation_date[i] = 0;
                legacy.defoliation_method[i] = 0;
                legacy.defoliant_app_rate[i] = 0.;
            }

            for ao in &self.agronomy_operations {
                match ao {
                    AgronomyOperation::irrigation(irrigation) => {
                        if irrigation.predict {
                            legacy.max_irrigation = irrigation.max_amount.unwrap();
                            legacy.day_start_pred_irrig = irrigation.date.ordinal() as i32;
                            legacy.day_stop_pred_irrig =
                                irrigation.stop_predict_date.unwrap().ordinal() as i32;
                            if let IrrigationMethod::Drip = irrigation.method {
                                legacy.location_column_drip = utils::slab_horizontal_location(
                                    irrigation.drip_x,
                                    legacy.row_space,
                                )?
                                    as i32;
                                legacy.location_layer_drip =
                                    utils::slab_vertical_location(irrigation.drip_y)? as i32;
                            }
                            legacy.irrig_method = irrigation.method as i32;
                        } else {
                            let idx = legacy.num_irrigations as usize;
                            let mut irrig = Irrig.write().expect("Irrig lock poisoned");
                            irrig[idx].day = irrigation.date.ordinal() as i32;
                            irrig[idx].amount = irrigation.amount;
                            if let IrrigationMethod::Drip = irrigation.method {
                                irrig[idx].LocationColumnDrip = utils::slab_horizontal_location(
                                    irrigation.drip_x,
                                    legacy.row_space,
                                )?
                                    as i32;
                                irrig[idx].LocationLayerDrip =
                                    utils::slab_vertical_location(irrigation.drip_y)? as i32;
                            }
                            irrig[idx].method = irrigation.method as i32;
                            legacy.num_irrigations += 1;
                        }
                    }
                    AgronomyOperation::defoliation {
                        date,
                        open_ratio,
                        predict,
                        ppa,
                    } => {
                        legacy.defoliation_date[idef] = date.ordinal() as i32;
                        legacy.defoliant_app_rate[idef] = if *predict { -99.9 } else { *ppa };
                        legacy.defoliation_method[idef] = *open_ratio;
                        legacy.day_first_def = legacy.defoliation_date[0];
                        idef += 1;
                    }
                    _ => {}
                }
            }
            self.num_watertable_data = self
                .agronomy_operations
                .iter()
                .filter(|x| matches!(x, AgronomyOperation::watertable { .. }))
                .count();
            if matches!(self.light_intercept_method, LightInterceptMethod::Latered)
                && self.light_intercept_parameters.is_none()
            {
                panic!(
                    "light_intercept_parameters must be provided when using the Latered light intercept method"
                );
            }
            init_soil(legacy, &self.soil_layers, &self.soil_hydraulic);
            initialize_root_data(legacy);
            //     initialize some variables at the start of simulation.
            legacy.soil_nitrogen_at_start =
                legacy.total_soil_no3_n + legacy.total_soil_nh4_n + legacy.total_soil_urea_n;
            legacy.plant_weight_at_start =
                legacy.total_root_weight + legacy.total_stem_weight + 0.2 + legacy.reserve_c;
            legacy.write_to_globals();
        }
        // If this is the first time the function is executed, get the ambient CO2 correction.
        self.ambient_CO2_factor = utils::ambient_CO2_factor(self.start_date.year());
        Ok(())
    }

    /// Write the output CSV file's header.
    pub fn output_file_headers(&self) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = csv::Writer::from_path(self.path.parent().unwrap().join("output.csv"))?;
        writer.write_field("date")?;
        writer.write_field("light_interception")?;
        writer.write_field("lint_yield")?;
        writer.write_field("leaf_area_index")?;
        writer.write_field("seed_cotton_yield")?;
        writer.write_field("plant_height")?;
        writer.write_field("main_stem_nodes")?;
        writer.write_field("leaf_weight")?;
        writer.write_field("petiole_weight")?;
        writer.write_field("stem_weight")?;
        writer.write_field("number_of_squares")?;
        writer.write_field("number_of_green_bolls")?;
        writer.write_field("number_of_open_bolls")?;
        writer.write_field("square_weight")?;
        writer.write_field("boll_weight")?;
        writer.write_field("root_weight")?;
        writer.write_field("plant_weight")?;
        writer.write_field("swc0-10")?;
        writer.write_field("swc0-20")?;
        writer.write_field("swc0-30")?;
        writer.write_field("swc1-10")?;
        writer.write_field("swc1-20")?;
        writer.write_field("swc1-30")?;
        writer.write_field("swc2-10")?;
        writer.write_field("swc2-20")?;
        writer.write_field("swc2-30")?;
        writer.write_field("swc3-10")?;
        writer.write_field("swc3-20")?;
        writer.write_field("swc3-30")?;
        writer.write_field("lai00")?;
        writer.write_field("lai01")?;
        writer.write_field("lai02")?;
        writer.write_field("lai03")?;
        writer.write_field("lai04")?;
        writer.write_field("lai05")?;
        writer.write_field("lai06")?;
        writer.write_field("lai07")?;
        writer.write_field("lai08")?;
        writer.write_field("lai09")?;
        writer.write_field("lai10")?;
        writer.write_field("lai11")?;
        writer.write_field("lai12")?;
        writer.write_field("lai13")?;
        writer.write_field("lai14")?;
        writer.write_field("lai15")?;
        writer.write_field("lai16")?;
        writer.write_field("lai17")?;
        writer.write_field("lai18")?;
        writer.write_field("lai19")?;
        writer.write_record(None::<&[u8]>)?;
        Ok(())
    }

    /// This function opens the soil root impedance data file and reads it. It is called from [Profile::initialize()]
    /// and executed once at the beginning of the simulation. The variables read here are later used to compute soil
    /// impedance to root growth.
    ///
    /// Returns the impedance curves that relate bulk density and soil water content to root
    /// penetration resistance. These tables are used by the root growth sub-model in Rust.
    ///
    fn read_soil_impedance(
        &self,
        path: &std::path::Path,
    ) -> Result<RootImpedanceTables, Box<dyn std::error::Error>> {
        let mut rdr = csv::Reader::from_path(path)?;
        let headers = rdr.headers()?.clone();
        let mut bulk_density = Vec::new();
        for value in headers.iter().skip(1) {
            bulk_density.push(
                value
                    .parse::<f64>()
                    .map_err(|e| Box::new(e) as Box<dyn std::error::Error>)?,
            );
        }
        let mut water_content = Vec::new();
        let mut impedance = vec![Vec::new(); bulk_density.len()];
        for result in rdr.records() {
            let record = result?;
            let water_content_value: f64 = record
                .get(0)
                .ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "missing water content column",
                    )
                })?
                .parse::<f64>()
                .map_err(|e| Box::new(e) as Box<dyn std::error::Error>)?;
            water_content.push(water_content_value);
            for (col_idx, item) in record.iter().enumerate().skip(1) {
                let value: f64 = item
                    .parse::<f64>()
                    .map_err(|e| Box::new(e) as Box<dyn std::error::Error>)?;
                impedance[col_idx - 1].push(value);
            }
        }
        Ok(RootImpedanceTables::new(
            water_content,
            bulk_density,
            impedance,
        ))
    }

    pub fn write_record(&self) -> Result<(), Box<dyn std::error::Error>> {
        let mut f = std::fs::OpenOptions::new()
            .write(true)
            .append(true)
            .open(self.path.parent().unwrap().join("output.csv"))?;

        let (
            date_string,
            light_intercept,
            lint_yield,
            leaf_area_index,
            seed_cotton_weight,
            plant_height,
            num_fruit_branches,
            total_leaf_weight,
            total_petiole_weight,
            total_stem_weight,
            num_squares,
            num_green_bolls,
            num_open_bolls,
            total_square_weight,
            boll_and_burr_weight,
            total_root_weight,
            above_ground_biomass,
            water_samples,
            leaf_area_indexes,
        ) = {
            let legacy = &self.model_state.legacy;
            let date_string = chrono::NaiveDate::from_yo_opt(legacy.iyear, legacy.daynum as u32)
                .unwrap()
                .format("%F")
                .to_string();
            let seed_cotton_weight = (legacy.cotton_weight_open_bolls
                + legacy.cotton_weight_green_bolls)
                * legacy.plant_population
                / 1000.;
            let boll_and_burr_weight = legacy.cotton_weight_open_bolls
                + legacy.cotton_weight_green_bolls
                + legacy.burr_weight_green_bolls
                + legacy.burr_weight_open_bolls;
            let above_ground_biomass = if legacy.daynum >= legacy.day_emerge && legacy.isw > 0 {
                (legacy.plant_weight - legacy.total_root_weight) * legacy.plant_population / 1000.
            } else {
                0.
            };
            let water_samples = [
                legacy.vol_water_content[[3, 0]],
                legacy.vol_water_content[[5, 0]],
                legacy.vol_water_content[[7, 0]],
                legacy.vol_water_content[[3, 4]],
                legacy.vol_water_content[[5, 4]],
                legacy.vol_water_content[[7, 4]],
                legacy.vol_water_content[[3, 8]],
                legacy.vol_water_content[[5, 8]],
                legacy.vol_water_content[[7, 8]],
                legacy.vol_water_content[[3, 12]],
                legacy.vol_water_content[[5, 12]],
                legacy.vol_water_content[[7, 12]],
            ];
            (
                date_string,
                legacy.light_intercept,
                legacy.lint_yield,
                legacy.leaf_area_index,
                seed_cotton_weight,
                legacy.plant_height,
                legacy.num_fruit_branches[0],
                self.model_state.total_leaf_weight(),
                legacy.total_petiole_weight,
                legacy.total_stem_weight,
                legacy.num_squares,
                legacy.num_green_bolls,
                legacy.num_open_bolls,
                legacy.total_square_weight,
                boll_and_burr_weight,
                legacy.total_root_weight,
                above_ground_biomass,
                water_samples,
                legacy.leaf_area_indexes.clone(),
            )
        };

        let mut record = vec![
            date_string,
            light_intercept.to_string(),
            lint_yield.to_string(),
            leaf_area_index.to_string(),
            seed_cotton_weight.to_string(),
            plant_height.to_string(),
            num_fruit_branches.to_string(),
            total_leaf_weight.to_string(),
            total_petiole_weight.to_string(),
            total_stem_weight.to_string(),
            num_squares.to_string(),
            num_green_bolls.to_string(),
            num_open_bolls.to_string(),
            total_square_weight.to_string(),
            boll_and_burr_weight.to_string(),
            total_root_weight.to_string(),
            above_ground_biomass.to_string(),
        ];
        for value in water_samples {
            record.push(value.to_string());
        }
        for value in leaf_area_indexes {
            record.push(value.to_string());
        }

        writeln!(f, "{}", record.join(","))?;
        Ok(())
    }
}
/// This function initializes many "global" variables at the start of a simulation.
///
/// It is called from [Profile::initialize()].
///
/// NOTE: that initialization is needed at the start of each simulation (NOT at start of the run).
fn initialize_global(legacy: &mut LegacyGlobalState) {
    legacy.abscised_fruit_sites = 0.;
    legacy.abscised_leaf_weight = 0.;
    legacy.addwtbl = 0.;
    legacy.applied_water = 0.;
    legacy.average_lwp = 0.;
    legacy.average_lwp_min = 0.;

    legacy.bloom_weight_loss = 0.;
    legacy.burr_n_conc = 0.;
    legacy.burr_nitrogen = 0.;
    legacy.burr_weight_green_bolls = 0.;
    legacy.burr_weight_open_bolls = 0.;

    legacy.carbon_allocated_for_root_growth = 0.;
    legacy.cotton_weight_green_bolls = 0.;
    legacy.cotton_weight_open_bolls = 0.;
    legacy.carbon_stress = 1.;
    legacy.cum_evaporation = 0.;
    legacy.cum_fertilizer_n = 0.;
    legacy.cum_net_photosynth = 0.;
    legacy.adj_square_absc = 0.;
    legacy.cum_nitrogen_uptake = 0.;
    legacy.cum_plant_n_loss = 0.;
    legacy.cum_transpiration = 0.;
    legacy.cum_water_added = 0.;
    legacy.cum_water_drained = 0.;

    legacy.extra_carbon = 0.;
    legacy.first_bloom = 0;
    legacy.first_square = 0;
    legacy.fruit_growth_ratio = 1.;

    legacy.ginp = 0.35;
    legacy.gintot = 0.35;
    legacy.green_bolls_lost = 0.;

    legacy.last_irrigation = 0;
    legacy.leaf_area_index = 0.001;
    legacy.leaf_n_conc = 0.056;
    legacy.leaf_nitrogen = 0.0112;
    legacy.lint_yield = 0.;

    legacy.max_irrigation = 0.;
    legacy.mineralized_organic_n = 0.;

    legacy.nitrogen_stress = 1.;
    legacy.num_abscised_leaves = 0;
    legacy.num_open_bolls = 0.;
    legacy.num_pre_fru_nodes = 1;
    legacy.num_shedding_tags = 0;
    legacy.num_veg_branches = 1;
    legacy.n_stress_fruiting = 1.;
    legacy.n_stress_roots = 1.;
    legacy.n_stress_veg = 1.;

    legacy.percent_defoliation = 0.;
    legacy.petiole_n_conc = 0.;
    legacy.petiole_nitrogen = 0.;
    legacy.petiole_no3_n_conc = 0.;
    legacy.pixcon = 0.;
    legacy.pixda = 1.;
    legacy.pixdn = 1.;
    legacy.pixdz = 1.;
    legacy.pix_in_plants = 0.;
    legacy.plant_height = 4.0;
    legacy.plant_weight = 0.;
    legacy.pot_gro_stem = 0.;

    legacy.reserve_c = 0.06;
    legacy.root_n_conc = 0.026;
    legacy.root_nitrogen = 0.0052;
    legacy.root_weight_loss = 0.;

    legacy.seed_n_conc = 0.;
    legacy.seed_nitrogen = 0.;
    legacy.soil_nitrogen_loss = 0.;
    legacy.square_n_conc = 0.;
    legacy.square_nitrogen = 0.;
    legacy.stem_n_conc = 0.036;
    legacy.stem_nitrogen = 0.0072;
    legacy.sum_no3_n90 = 0.;
    legacy.supply_nh4_n = 0.;
    legacy.supply_no3_n = 0.;

    legacy.total_actual_leaf_growth = 0.;
    legacy.total_actual_petiole_growth = 0.;
    legacy.total_petiole_weight = 0.;
    legacy.total_required_n = 0.;
    legacy.total_square_weight = 0.;
    legacy.total_stem_weight = 0.2;

    legacy.water_stress = 1.;
    legacy.water_stress_stem = 1.;
    legacy.water_table_layer = 1000;
    //
    for i in 0..3 {
        legacy.delay_new_fru_branch[i] = 0.;
        legacy.lwp_min_x[i] = 0.;
        legacy.lwp_x[i] = 0.;
        legacy.num_fruit_branches[i] = 0;
        for j in 0..30 {
            legacy.delay_new_node[[i, j]] = 0.;
            legacy.leaf_area_main_stem[[i, j]] = 0.;
            legacy.leaf_weight_main_stem[[i, j]] = 0.;
            legacy.num_nodes[[i, j]] = 0;
            legacy.petiole_weight_main_stem[[i, j]] = 0.;
            legacy.pot_gro_leaf_area_main_stem[[i, j]] = 0.;
            legacy.pot_gro_leaf_weight_main_stem[[i, j]] = 0.;
            legacy.pot_gro_petiole_weight_main_stem[[i, j]] = 0.;
            legacy.node_layer[[i, j]] = 0;
            for k in 0..5 {
                legacy.age_of_boll[[i, j, k]] = 0.;
                legacy.age_of_site[[i, j, k]] = 0.;
                legacy.avrg_node_temper[[i, j, k]] = 0.;
                legacy.boll_weight[[i, j, k]] = 0.;
                legacy.burr_weight[[i, j, k]] = 0.;
                legacy.fruiting_code[[i, j, k]] = 0;
                legacy.fruit_fraction[[i, j, k]] = 0.;
                legacy.leaf_age[[i, j, k]] = 0.;
                legacy.leaf_area_nodes[[i, j, k]] = 0.;
                legacy.leaf_weight_nodes[[i, j, k]] = 0.;
                legacy.petiole_weight_nodes[[i, j, k]] = 0.;
                legacy.pot_gro_bolls[[i, j, k]] = 0.;
                legacy.pot_gro_burrs[[i, j, k]] = 0.;
                legacy.pot_gro_leaf_area_nodes[[i, j, k]] = 0.;
                legacy.pot_gro_leaf_weight_nodes[[i, j, k]] = 0.;
                legacy.pot_gro_petiole_weight_nodes[[i, j, k]] = 0.;
                legacy.pot_gro_squares[[i, j, k]] = 0.;
                legacy.square_weight[[i, j, k]] = 0.;
            }
        }
    }
    //
    for i in 0..9 {
        legacy.age_of_pre_fru_node[i] = 0.;
        legacy.leaf_area_pre_fru[i] = 0.;
        legacy.pot_gro_leaf_area_pre_fru[i] = 0.;
        legacy.pot_gro_leaf_weight_pre_fru[i] = 0.;
        legacy.pot_gro_petiole_weight_pre_fru[i] = 0.;
        legacy.leaf_weight_pre_fru[i] = 0.;
        legacy.petiole_weight_pre_fru[i] = 0.;
        legacy.node_layer_pre_fru[i] = 0;
    }
    //
    for i in 0..20 {
        legacy.abscission_lag[i] = 0.;
        legacy.shed_by_carbon_stress[i] = 0.;
        legacy.shed_by_nitrogen_stress[i] = 0.;
        legacy.shed_by_water_stress[i] = 0.;
    }
    for i in 0..20 {
        legacy.leaf_area[i] = 0.;
        legacy.leaf_weight_layer[i] = 0.;
    }
    legacy.leaf_weight_layer[0] = 0.2;
    //
    for k in 0..maxk as usize {
        legacy.foliage_temp[k] = 295.;
        legacy.mulch_temp[k] = 295.;
    }
    //
    for l in 0..maxl as usize {
        legacy.rlat1[l] = 0.;
        legacy.rlat2[l] = 0.;
        for k in 0..maxk as usize {
            legacy.root_wt_capbl_uptake[[l, k]] = 0.;
            legacy.root_impede[[l, k]] = 0.;
        }
    }
    //
    for i in 0..365 {
        legacy.stem_weight[i] = 0.;
    }
}

/// This function initializes the root submodel parameters and variables.
/// It is called by ReadInput(). it is executed once at the beginning of the simulation.
///
/// Global or file scope variables referenced:
/// dl, PlantRowColumn, nk, nl, PerPlantArea, RowSpace.
///
/// Global or file scope variables set:
/// ActualRootGrowth[maxl][maxk], cgind[3], DepthLastRootLayer,
/// LastTaprootLayer, LateralRootFlag[maxl], NumLayersWithRoots, NumRootAgeGroups,
/// PotGroRoots[maxl][maxk], RootAge[maxl][maxk], RootColNumLeft[maxl],
/// RootColNumRight[maxl], RootGroFactor[maxl][maxk],
/// RootWeight[maxl][maxk][3], TapRootLength, TotalRootWeight.
fn initialize_root_data(legacy: &mut LegacyGlobalState) {
    // The parameters of the root model are defined for each root class:
    // grind(i), cuind(i), thtrn(i), trn(i), thdth(i), dth(i).
    legacy.num_root_age_groups = 3;
    legacy.cgind[0] = 1.;
    legacy.cgind[1] = 1.;
    legacy.cgind[2] = 0.10;
    let rlint = 10.; // Vertical interval, in cm, along the taproot, for initiating lateral roots.
    let mut ll = 1; // Counter for layers with lateral roots.
    let mut sumdl = 0.; // Distance from soil surface to the middle of a soil layer.
    for l in 0..legacy.nl as usize {
        // Using the value of rlint (interval between lateral roots), the layers from which lateral roots may be initiated are now computed.
        // legacy.lateral_root_flag[l] is assigned a value of 1 for these layers.
        legacy.lateral_root_flag[l] = 0;
        if l > 0 {
            sumdl += 0.5 * legacy.dl[l - 1];
        }
        sumdl += 0.5 * legacy.dl[l];
        if sumdl >= ll as f64 * rlint {
            legacy.lateral_root_flag[l] = 1;
            ll += 1;
        }
    }
    // All the state variables of the root system are initialized to zero.
    for l in 0..legacy.nl as usize {
        if l < 3 {
            legacy.root_col_num_left[l] = legacy.plant_row_column - 1;
            legacy.root_col_num_right[l] = legacy.plant_row_column + 2;
        } else if l < 7 {
            legacy.root_col_num_left[l] = legacy.plant_row_column;
            legacy.root_col_num_right[l] = legacy.plant_row_column + 1;
        } else {
            legacy.root_col_num_left[l] = 0;
            legacy.root_col_num_right[l] = 0;
        }
        //
        for k in 0..legacy.nk as usize {
            legacy.pot_gro_roots[[l, k]] = 0.;
            legacy.root_gro_factor[[l, k]] = 1.;
            legacy.actual_root_growth[[l, k]] = 0.;
            legacy.root_age[[l, k]] = 0.;
            for i in 0..3 {
                legacy.root_weight[[l, k, i]] = 0.;
            }
        }
    }
    //
    legacy.root_weight[[0, (legacy.plant_row_column - 1) as usize, 0]] = 0.0020;
    legacy.root_weight[[0, legacy.plant_row_column as usize, 0]] = 0.0070;
    legacy.root_weight[[0, legacy.plant_row_column as usize + 1, 0]] = 0.0070;
    legacy.root_weight[[0, legacy.plant_row_column as usize + 2, 0]] = 0.0020;
    legacy.root_weight[[1, (legacy.plant_row_column - 1) as usize, 0]] = 0.0040;
    legacy.root_weight[[1, legacy.plant_row_column as usize, 0]] = 0.0140;
    legacy.root_weight[[1, legacy.plant_row_column as usize + 1, 0]] = 0.0140;
    legacy.root_weight[[1, legacy.plant_row_column as usize + 2, 0]] = 0.0040;
    legacy.root_weight[[2, (legacy.plant_row_column - 1) as usize, 0]] = 0.0060;
    legacy.root_weight[[2, legacy.plant_row_column as usize, 0]] = 0.0210;
    legacy.root_weight[[2, legacy.plant_row_column as usize + 1, 0]] = 0.0210;
    legacy.root_weight[[2, legacy.plant_row_column as usize + 2, 0]] = 0.0060;
    legacy.root_weight[[3, legacy.plant_row_column as usize, 0]] = 0.0200;
    legacy.root_weight[[3, legacy.plant_row_column as usize + 1, 0]] = 0.0200;
    legacy.root_weight[[4, legacy.plant_row_column as usize, 0]] = 0.0150;
    legacy.root_weight[[4, legacy.plant_row_column as usize + 1, 0]] = 0.0150;
    legacy.root_weight[[5, legacy.plant_row_column as usize, 0]] = 0.0100;
    legacy.root_weight[[5, legacy.plant_row_column as usize + 1, 0]] = 0.0100;
    legacy.root_weight[[6, legacy.plant_row_column as usize, 0]] = 0.0050;
    legacy.root_weight[[6, legacy.plant_row_column as usize + 1, 0]] = 0.0050;
    // Start loop for all soil layers containing roots.
    legacy.depth_last_root_layer = 0.;
    legacy.total_root_weight = 0.;
    for l in 0..7 {
        legacy.depth_last_root_layer += legacy.dl[l]; // compute total depth to the last layer with roots (legacy.depth_last_root_layer).
        for k in 0..legacy.nk as usize {
            // For each soil soil cell with roots, compute total root weight per plant (legacy.total_root_weight), and convert legacy.root_weight from g per plant to g per cell.
            for i in 0..3 {
                legacy.total_root_weight += legacy.root_weight[[l, k, i]];
                legacy.root_weight[[l, k, i]] =
                    legacy.root_weight[[l, k, i]] * 0.01 * legacy.row_space / legacy.per_plant_area;
            }
            // initialize legacy.root_age to a non-zero value for each cell containing roots.
            if legacy.root_weight[[l, k, 0]] > 0. {
                legacy.root_age[[l, k]] = 0.01;
            }
        }
    }
    //     Initial value of taproot length, legacy.tap_root_length, is computed to the
    // middle of the last layer with roots. The last soil layer with
    // taproot, legacy.last_taproot_layer, is defined.
    legacy.num_layers_with_roots = 7;
    legacy.tap_root_length =
        legacy.depth_last_root_layer - 0.5 * legacy.dl[(legacy.num_layers_with_roots - 1) as usize];
    legacy.last_taproot_layer = 6;
}
