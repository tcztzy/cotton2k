//! # Cotton2K
//!
//! Simulation model for cotton

#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(dead_code)]
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
mod adjust;
mod utils;
use chrono::Datelike;
use chrono::NaiveDate;
use serde::Deserialize;
use std::io::Write;
use std::path::PathBuf;

#[derive(Debug)]
pub struct Cotton2KError {
    level: u8,
    message: String,
}

impl std::fmt::Display for Cotton2KError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", "Cotton2KError")
    }
}

impl std::error::Error for Cotton2KError {}

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

#[cfg(target_os = "windows")]
#[inline]
fn zero_i32() -> i32 {
    0
}

#[cfg(target_os = "windows")]
#[derive(Deserialize, Debug)]
pub struct Profile {
    #[serde(skip)]
    pub path: PathBuf,
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
    /// day after emergence of the previous plant map adjustment.
    #[serde(skip)]
    kprevadj: u32,
    /// maximum leaf area index.
    #[serde(skip)]
    lmax: f64,
    /// The effect of moisture stress on the photosynthetic rate
    #[serde(skip)]
    ptsred: f64,
    /// correction factor for ambient CO2 in air
    #[serde(skip)]
    ambient_CO2_factor: f64,
    #[serde(skip)]
    states: Vec<State>,
}

#[cfg(target_os = "linux")]
#[inline]
fn zero_u32() -> u32 {
    0
}

#[cfg(target_os = "linux")]
#[derive(Deserialize, Debug)]
pub struct Profile {
    #[serde(skip)]
    pub path: PathBuf,
    #[serde(default = "zero_u32")]
    pub light_intercept_method: u32,
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
    /// day after emergence of the previous plant map adjustment.
    #[serde(skip)]
    kprevadj: u32,
    /// maximum leaf area index.
    #[serde(skip)]
    lmax: f64,
    /// The effect of moisture stress on the photosynthetic rate
    #[serde(skip)]
    ptsred: f64,
    /// correction factor for ambient CO2 in air
    #[serde(skip)]
    ambient_CO2_factor: f64,
    #[serde(skip)]
    states: Vec<State>,
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
            start_date: NaiveDate::from_ymd(1900, 1, 1),
            stop_date: NaiveDate::from_ymd(2100, 1, 1),
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
pub enum PixMethod {
    Banded = 0,
    Sprinkler = 1,
    Broadcast = 2,
}

#[derive(Deserialize, Debug, Clone, Copy)]
#[serde(tag = "type")]
pub enum AgronomyOperation {
    irrigation {
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
        method: PixMethod,
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

#[derive(Debug)]
pub struct State {
    pub plant_height: f64,
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

/// This function initializes many "global" variables at the start of a simulation.
///
/// It is called from [Profile::initialize()].
///
/// NOTE: that initialization is needed at the start of each simulation (NOT at start of the run).
unsafe fn InitializeGlobal() {
    AbscisedFruitSites = 0.;
    AbscisedLeafWeight = 0.;
    addwtbl = 0.;
    AppliedWater = 0.;
    AverageLwp = 0.;
    AverageLwpMin = 0.;

    BloomWeightLoss = 0.;
    BurrNConc = 0.;
    BurrNitrogen = 0.;
    BurrWeightGreenBolls = 0.;
    BurrWeightOpenBolls = 0.;

    CarbonAllocatedForRootGrowth = 0.;
    CottonWeightGreenBolls = 0.;
    CottonWeightOpenBolls = 0.;
    CarbonStress = 1.;
    CumEvaporation = 0.;
    CumFertilizerN = 0.;
    CumNetPhotosynth = 0.;
    AdjSquareAbsc = 0.;
    CumNitrogenUptake = 0.;
    CumPlantNLoss = 0.;
    CumTranspiration = 0.;
    CumWaterAdded = 0.;
    CumWaterDrained = 0.;

    ExtraCarbon = 0.;
    FirstBloom = 0;
    FirstSquare = 0;
    FruitGrowthRatio = 1.;

    ginp = 0.35;
    Gintot = 0.35;
    GreenBollsLost = 0.;

    LastIrrigation = 0;
    LeafAreaIndex = 0.001;
    LeafNConc = 0.056;
    LeafNitrogen = 0.0112;
    LintYield = 0.;

    MaxIrrigation = 0.;
    MineralizedOrganicN = 0.;

    NitrogenStress = 1.;
    NumAbscisedLeaves = 0;
    NumOpenBolls = 0.;
    NumPreFruNodes = 1;
    NumSheddingTags = 0;
    NumVegBranches = 1;
    NumWaterTableData = 0;
    NStressFruiting = 1.;
    NStressRoots = 1.;
    NStressVeg = 1.;

    PercentDefoliation = 0.;
    PetioleNConc = 0.;
    PetioleNitrogen = 0.;
    PetioleNO3NConc = 0.;
    pixcon = 0.;
    pixda = 1.;
    pixdn = 1.;
    pixdz = 1.;
    PixInPlants = 0.;
    PlantHeight = 4.0;
    PlantWeight = 0.;
    PotGroStem = 0.;

    ReserveC = 0.06;
    RootNConc = 0.026;
    RootNitrogen = 0.0052;
    RootWeightLoss = 0.;

    SeedNConc = 0.;
    SeedNitrogen = 0.;
    SoilNitrogenLoss = 0.;
    SquareNConc = 0.;
    SquareNitrogen = 0.;
    StemNConc = 0.036;
    StemNitrogen = 0.0072;
    SumNO3N90 = 0.;
    SupplyNH4N = 0.;
    SupplyNO3N = 0.;

    TotalActualLeafGrowth = 0.;
    TotalActualPetioleGrowth = 0.;
    TotalPetioleWeight = 0.;
    TotalRequiredN = 0.;
    TotalSquareWeight = 0.;
    TotalStemWeight = 0.2;

    WaterStress = 1.;
    WaterStressStem = 1.;
    WaterTableLayer = 1000;
    //
    for i in 0..3 {
        DelayNewFruBranch[i] = 0.;
        LwpMinX[i] = 0.;
        LwpX[i] = 0.;
        NumFruitBranches[i] = 0;
        for j in 0..30 {
            DelayNewNode[i][j] = 0.;
            LeafAreaMainStem[i][j] = 0.;
            LeafWeightMainStem[i][j] = 0.;
            NumNodes[i][j] = 0;
            PetioleWeightMainStem[i][j] = 0.;
            PotGroLeafAreaMainStem[i][j] = 0.;
            PotGroLeafWeightMainStem[i][j] = 0.;
            PotGroPetioleWeightMainStem[i][j] = 0.;
            NodeLayer[i][j] = 0;
            for k in 0..5 {
                AgeOfBoll[i][j][k] = 0.;
                AgeOfSite[i][j][k] = 0.;
                AvrgNodeTemper[i][j][k] = 0.;
                BollWeight[i][j][k] = 0.;
                BurrWeight[i][j][k] = 0.;
                FruitingCode[i][j][k] = 0;
                FruitFraction[i][j][k] = 0.;
                LeafAge[i][j][k] = 0.;
                LeafAreaNodes[i][j][k] = 0.;
                LeafWeightNodes[i][j][k] = 0.;
                PetioleWeightNodes[i][j][k] = 0.;
                PotGroBolls[i][j][k] = 0.;
                PotGroBurrs[i][j][k] = 0.;
                PotGroLeafAreaNodes[i][j][k] = 0.;
                PotGroLeafWeightNodes[i][j][k] = 0.;
                PotGroPetioleWeightNodes[i][j][k] = 0.;
                PotGroSquares[i][j][k] = 0.;
                SquareWeight[i][j][k] = 0.;
            }
        }
    }
    //
    for i in 0..9 {
        AgeOfPreFruNode[i] = 0.;
        LeafAreaPreFru[i] = 0.;
        PotGroLeafAreaPreFru[i] = 0.;
        PotGroLeafWeightPreFru[i] = 0.;
        PotGroPetioleWeightPreFru[i] = 0.;
        LeafWeightPreFru[i] = 0.;
        PetioleWeightPreFru[i] = 0.;
        NodeLayerPreFru[i] = 0;
    }
    //
    for i in 0..20 {
        AbscissionLag[i] = 0.;
        DayWaterTableInput[i] = 0;
        ElCondSatSoil[i] = 0.;
        LevelsOfWaterTable[i] = 0.;
        ShedByCarbonStress[i] = 0.;
        ShedByNitrogenStress[i] = 0.;
        ShedByWaterStress[i] = 0.;
    }
    for i in 0..20 {
        LeafArea[i] = 0.;
        LeafWeightLayer[i] = 0.;
    }
    LeafWeightLayer[0] = 0.2;
    //
    for k in 0..maxk as usize {
        FoliageTemp[k] = 295.;
        MulchTemp[k] = 295.;
    }
    //
    for l in 0..maxl as usize {
        rlat1[l] = 0.;
        rlat2[l] = 0.;
        for k in 0..maxk as usize {
            RootWtCapblUptake[l][k] = 0.;
            RootImpede[l][k] = 0.;
        }
    }
    //
    for i in 0..365 {
        StemWeight[i] = 0.;
    }
}

/// This function opens the initial soil data file and reads it.
/// It is executed once at the beginning of the simulation.
/// It is called by [Profile::initialize()].
unsafe fn InitSoil(soil_layers: &[SoilLayer; 14], soil_hydraulic: &SoilHydraulic) {
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
        airdr[i] = 0.;
        thetas[i] = 0.;
        alpha[i] = 0.;
        beta[i] = 0.;
        SaturatedHydCond[i] = 0.;
        BulkDensity[i] = 0.;
    }
    RatioImplicit = soil_hydraulic.implicit_ratio;
    conmax = soil_hydraulic.max_conductivity;
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
        airdr[il] = layer.theta_d;
        thetas[il] = layer.theta_s;
        alpha[il] = layer.alpha;
        beta[il] = layer.beta;
        SaturatedHydCond[il] = layer.hcs;
        BulkDensity[il] = layer.bulk_density;
    }

    let mut j = 0usize; // horizon number
    let mut sumdl = 0f64; // depth to the bottom this layer (cm);
    let rm = 2.65f64; // density of the solid fraction of the soil (g / cm3)
    let mut bdl = [0f64; 40]; // array of bulk density of soil layers

    for l in 0..nl as usize {
        //     Using the depth of each horizon layer (ldepth), the horizon
        //  number (SoilHorizonNum) is computed for each soil layer.
        sumdl += dl[l];
        while (sumdl > ldepth[j]) && (j < lyrsol) {
            j += 1;
        }
        SoilHorizonNum[l] = j as i32;
        // bdl, thad, thts are defined for each soil layer, using the respective input variables BulkDensity, airdr,
        // thetas.
        //
        // FieldCapacity, MaxWaterCapacity and thetar are computed for each layer, as water content ($cm^3\ cm^{-3}$)
        // of each layer corresponding to matric potentials of psisfc (for field capacity), psidra (for free drainage)
        // and -15 bars (for permanent wilting point), respectively, using function qpsi.
        //
        // pore space volume (PoreSpace) is also computed for each layer. make sure that saturated water content is not
        // more than pore space.
        bdl[l] = BulkDensity[j];
        PoreSpace[l] = 1. - BulkDensity[j] / rm;
        if thetas[j] > PoreSpace[l] {
            thetas[j] = PoreSpace[l];
        }
        thad[l] = airdr[j];
        thts[l] = thetas[j];
        FieldCapacity[l] = qpsi(psisfc, thad[l], thts[l], alpha[j], beta[j]);
        MaxWaterCapacity[l] = qpsi(psidra, thad[l], thts[l], alpha[j], beta[j]);
        thetar[l] = qpsi(-15., thad[l], thts[l], alpha[j], beta[j]);
        // When the saturated hydraulic conductivity (SaturatedHydCond) is not given, it is computed from the hydraulic
        // conductivity at field capacity (condfc), using the wcond function.
        if SaturatedHydCond[j] <= 0. {
            SaturatedHydCond[j] =
                condfc[j] / wcond(FieldCapacity[l], thad[l], thts[l], beta[j], 1., 1.);
        }
    }
    // Loop for all soil layers. Compute depth from soil surface to the end of each layer (sumdl).
    sumdl = 0.;
    for l in 0..nl as usize {
        sumdl += dl[l];
        // At start of simulation compute estimated movable fraction of nitrates in each soil layer, following the work
        // of:
        //     Bowen, W.T., Jones, J.W., Carsky, R.J., and Quintana, J.O. 1993.
        //  Evaluation of the nitrogen submodel of CERES-maize following legume
        //  green manure incorporation. Agron. J. 85:153-159.
        //
        // The fraction of total nitrate in a layer that is in solution and can move from one layer to the next with
        // the downward flow of water, FLOWNO3[l], is a function of the adsorption coefficient, soil bulk density, and
        // the volumetric soil water content at the drained upper limit.
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
        NO3FlowFraction[l] = 1. / (1. + coeff * bdl[l] / MaxWaterCapacity[l]);
        // Determine the corresponding 15 cm layer of the input file. Compute the initial volumetric water content
        // (VolWaterContent) of each layer, and check that it will not be less than the air-dry value or more than pore
        // space volume.
        j = ((sumdl - 1.) / 15.).floor() as usize;
        if j > 13 {
            j = 13;
        }
        let n = SoilHorizonNum[l] as usize;
        VolWaterContent[l][0] = FieldCapacity[l] * h2oint[j] / 100.;
        if VolWaterContent[l][0] < airdr[n] {
            VolWaterContent[l][0] = airdr[n];
        }
        if VolWaterContent[l][0] > PoreSpace[l] {
            VolWaterContent[l][0] = PoreSpace[l];
        }
        // Initial values of ammonium N (rnnh4, VolNh4NContent) and nitrate N (rnno3, VolNo3NContent) are converted
        // from kgs per ha to $mg\ cm^{-3}$ for each soil layer, after checking for minimal amounts.
        if rnno3[j] < 2.0 {
            rnno3[j] = 2.0;
        }
        if rnnh4[j] < 0.2 {
            rnnh4[j] = 0.2;
        }
        VolNo3NContent[l][0] = rnno3[j] / 15. * 0.01;
        VolNh4NContent[l][0] = rnnh4[j] / 15. * 0.01;
        // organic matter in mg / cm3 units.
        let om = (oma[j] / 100.) * bdl[l] * 1000.;
        // potom is the proportion of readily mineralizable om. it is a function of soil depth (sumdl, in cm), modified
        // from GOSSYM (where it probably includes the 0.4 factor for organic C in om).
        let mut potom = 0.15125 - 0.02878 * sumdl.ln();
        if potom < 0. {
            potom = 0.;
        }
        // FreshOrganicMatter is the readily mineralizable organic matter (="fresh organic matter" in CERES models).
        // HumusOrganicMatter is the remaining organic matter, which is mineralized very slowly.
        FreshOrganicMatter[l][0] = om * potom;
        HumusOrganicMatter[l][0] = om * (1. - potom);
    }
    // Since the initial value has been set for the first column only in each layer, these values are now assigned to
    // all the other columns.
    for l in 0..nl as usize {
        for k in 1..nk as usize {
            VolWaterContent[l][k] = VolWaterContent[l][0];
            VolNo3NContent[l][k] = VolNo3NContent[l][0];
            VolNh4NContent[l][k] = VolNh4NContent[l][0];
            FreshOrganicMatter[l][k] = FreshOrganicMatter[l][0];
            HumusOrganicMatter[l][k] = HumusOrganicMatter[l][0];
        }
    }
    // Total amounts of water (InitialTotalSoilWater), nitrate N (TotalSoilNo3N), ammonium N (TotalSoilNh4N), and urea
    // N (TotalSoilUreaN) are computed for the whole slab.
    InitialTotalSoilWater = 0.;
    TotalSoilNo3N = 0.;
    TotalSoilNh4N = 0.;
    TotalSoilUreaN = 0.;

    for l in 0..nl as usize {
        for k in 0..nk as usize {
            InitialTotalSoilWater += VolWaterContent[l][k] * dl[l] * wk[k];
            TotalSoilNo3N += VolNo3NContent[l][k] * dl[l] * wk[k];
            TotalSoilNh4N += VolNh4NContent[l][k] * dl[l] * wk[k];
            VolUreaNContent[l][k] = 0.;
        }
    }
    // InitialTotalSoilWater is converted from cm3 per slab to mm.
    InitialTotalSoilWater = 10. * InitialTotalSoilWater / RowSpace;
    let bsand = 20f64; // heat conductivity of sand and silt (mcal cm-1 s-1 C-1).
    let bclay = 7f64; // heat conductivity of clay (mcal cm-1 s-1 C-1).
    let cka = 0.0615f64; // heat conductivity of air (mcal cm-1 s-1 C-1).
    let ckw = 1.45f64; // heat conductivity of water (mcal cm-1 s-1 C-1).
    let cmin = 0.46f64; // heat capacity of the mineral fraction of the soil.
    let corg = 0.6f64; // heat capacity of the organic fraction of the soil.
    let ga = 0.144f64; // shape factor for air in pore spaces.
    let ro = 1.3f64; // specific weight of organic fraction of soil.

    // Compute aggregation factors:
    dsand = form(bsand, ckw, ga); // aggregation factor for sand in water
    dclay = form(bclay, ckw, ga); // aggregation factor for clay in water
    let dsandair: f64 = form(bsand, cka, ga); // aggregation factor for sand in air
    let dclayair: f64 = form(bclay, cka, ga); // aggregation factor for clay in air

    // Loop over all soil layers, and define indices for some soil arrays.

    sumdl = 0.; // sum of depth of consecutive soil layers.

    for l in 0..nl as usize {
        sumdl += dl[l];
        let mut j = ((sumdl + 14.) / 15.).floor() as usize - 1; // layer definition for oma
        if j > 13 {
            j = 13;
        }
        // Using the values of the clay and organic matter percentages in the soil, compute mineral and organic
        // fractions of the soil, by weight and by volume.
        let mmo = oma[j] / 100.; // organic matter fraction of dry soil (by weight).
        let mm = 1. - mmo; // mineral fraction of dry soil (by weight).

        // MarginalWaterContent is set as a function of the sand fraction of the soil.
        let ra = (mmo / ro) / (mm / rm); // volume ratio of organic to mineral soil fractions.

        let i1 = SoilHorizonNum[l] as usize; //  layer definition as in soil hydrology input file.

        // The volume fractions of clay (ClayVolumeFraction) and of sand plus silt (SandVolumeFraction), are calculated
        MarginalWaterContent[l] = 0.1 - 0.07 * psand[i1] / 100.;
        let xo = (1. - PoreSpace[l]) * ra / (1. + ra); // organic fraction of soil (by volume).
        let xm = (1. - PoreSpace[l]) - xo; // mineral fraction of soil (by volume).
        ClayVolumeFraction[l] = pclay[i1] * xm / mm / 100.;
        SandVolumeFraction[l] = 1. - PoreSpace[l] - ClayVolumeFraction[l];
        // Heat capacity of the solid soil fractions (mineral + organic, by volume )
        HeatCapacitySoilSolid[l] = xm * cmin + xo * corg;
        // The heat conductivity of dry soil (HeatCondDrySoil) is computed using the procedure suggested by De Vries.
        HeatCondDrySoil[l] = 1.25
            * (PoreSpace[l] * cka
                + dsandair * bsand * SandVolumeFraction[l]
                + dclayair * bclay * ClayVolumeFraction[l])
            / (PoreSpace[l] + dsandair * SandVolumeFraction[l] + dclayair * ClayVolumeFraction[l]);
    }
}

fn drop_leaf_age(lai: f64) -> f64 {
    140. - 1. * lai
}

impl Profile {
    /// Run this profile.
    pub fn run(self: &mut Self) -> Result<(), Box<dyn std::error::Error>> {
        self.initialize()?;
        self.output_file_headers()?;
        unsafe {
            // Do daily simulations
            Daynum = DayStart - 1;
            bEnd = false;
            // Start the daily loop. If variable bEnd has been assigned a value of true end simulation.
            for _ in DayStart..(DayFinish + 1) {
                let bAdjustToDo = self.adjust()?;
                // Execute simulation for this day.
                match self.simulate_this_day() {
                    Err(e) => {
                        if e.level == 0 {
                            println!("{}", e.message);
                            break;
                        }
                    }
                    _ => {}
                }
                self.write_record()?;
                // If there are pending plant adjustments, call WriteStateVariables() to write
                // state variables of this day in a scratch file.
                if bAdjustToDo {
                    WriteStateVariables(true);
                }
                if bEnd {
                    break;
                }
            }
        }
        Ok(())
    }

    pub fn initialize(self: &mut Self) -> Result<(), Box<dyn std::error::Error>> {
        self.kprevadj = 0;
        unsafe {
            InitializeGlobal();
            light_intercept_method = self.light_intercept_method;
            Latitude = self.latitude;
            Longitude = self.longitude;
            Elevation = self.elevation;
            iyear = self.start_date.year();
            DayStart = self.start_date.ordinal() as i32;
            DayFinish = self.stop_date.ordinal() as i32;
            match self.emerge_date {
                Some(date) => DayEmerge = date.ordinal() as i32,
                None => {}
            }
            match self.plant_date {
                Some(date) => {
                    DayPlant = date.ordinal() as i32;
                }
                None => {}
            }
            if self.emerge_date.is_none() {
                // If the date of emergence has not been given, emergence will be simulated by the model. In this case,
                // isw = 0, and a check is performed to make sure that the date of planting has been given.
                if self.plant_date.is_none() {
                    panic!(
                        "one of planting date or emergence date must be given in the profile file!!"
                    );
                }
                isw = 0;
            } else if self.emerge_date > self.plant_date {
                // If the date of emergence has been given in the input: isw = 1 if simulation starts before emergence,
                isw = 1;
            } else {
                // or isw = 2 if simulation starts at emergence.
                isw = 2;
                Kday = 1;
            }
            // For advanced users only: if there is CO2 enrichment, read also CO2 factor, DOY dates
            // If soil mulch is used, read relevant parameters.
            match self.mulch {
                Some(mulch) => match mulch.indicator {
                    MulchType::NoMulch => {}
                    _ => {
                        MulchIndicator = mulch.indicator as i32;
                        MulchTranSW = mulch.sw_trans;
                        MulchTranLW = mulch.lw_trans;
                        DayStartMulch = mulch.start_date.ordinal() as i32;
                        match mulch.stop_date {
                            Some(date) => DayEndMulch = date.ordinal() as i32,
                            None => DayEndMulch = DayFinish,
                        }
                    }
                },
                None => {
                    MulchIndicator = MulchType::NoMulch as i32;
                }
            }
            SitePar[1] = self.site.wind_blow_after_sunrise;
            SitePar[2] = self.site.wind_max_after_noon;
            SitePar[3] = self.site.wind_stop_after_sunset;
            SitePar[4] = self.site.night_time_wind_factor;
            SitePar[7] = self.site.cloud_type_correction_factor;
            SitePar[8] = self.site.max_temperature_after_noon;
            SitePar[9] = self.site.deep_soil_temperature.0;
            SitePar[10] = self.site.deep_soil_temperature.1;
            SitePar[11] = self.site.deep_soil_temperature.2;
            SitePar[12] = self.site.dew_point_range.0;
            SitePar[13] = self.site.dew_point_range.1;
            SitePar[14] = self.site.dew_point_range.2;
            SitePar[15] = self.site.albedo_range.1;
            SitePar[16] = self.site.albedo_range.0;
            for pair in self.cultivar_parameters.iter().enumerate() {
                VarPar[pair.0 + 1] = *pair.1;
            }
            RowSpace = if self.skip_row_width > 0. {
                (self.row_space + self.skip_row_width) / 2.
            } else {
                self.row_space
            };
            // PlantRowLocation is the distance from edge of slab, cm, of the plant row.
            PlantRowLocation = RowSpace / 2.;
            // Compute PlantPopulation - number of plants per hectar, and PerPlantArea - the average surface area per
            // plant, in $dm^2$, and the empirical plant density factor (DensityFactor). This factor will be used to
            // express the effect of plant density on some plant growth rate functions.
            // Note that DensityFactor =1 for 5 plants per sq m (or 50000 per ha).
            PlantPopulation = self.plants_per_meter / RowSpace * 1000000.;
            PerPlantArea = 1000000. / PlantPopulation;
            DensityFactor = (VarPar[1] * (5. - PlantPopulation / 10000.)).exp();
            // Define the numbers of rows and columns in the soil slab (nl, nk).
            // Define the depth, in cm, of consecutive nl layers.
            nl = maxl;
            nk = maxk;
            dl[0] = 2.;
            dl[1] = 2.;
            dl[2] = 2.;
            dl[3] = 4.;
            for i in 4..(maxl - 2) as usize {
                dl[i] = 5.;
            }
            dl[(maxl - 2) as usize] = 10.;
            dl[(maxl - 1) as usize] = 10.;
            //      The width of the slab columns is computed by dividing the row
            //  spacing by the number of columns. It is assumed that slab width is
            //  equal to the average row spacing, and column widths are uniform.
            //      Note: wk is an array - to enable the option of non-uniform
            //  column widths in the future.
            //      PlantRowColumn (the column including the plant row) is now computed
            //      from
            //  PlantRowLocation (the distance of the plant row from the edge of the
            //  slab).
            let mut sumwk = 0.; // sum of column widths
            PlantRowColumn = 0;
            for k in 0..nk {
                wk[k as usize] = RowSpace / nk as f64;
                sumwk = sumwk + wk[k as usize];
                if PlantRowColumn == 0 && sumwk > PlantRowLocation {
                    PlantRowColumn = if (sumwk - PlantRowLocation) > (0.5 * wk[k as usize]) {
                        k - 1
                    } else {
                        k
                    };
                }
            }
        }
        let mut rdr = csv::Reader::from_path(&self.weather_path)?;
        let mut jdd: i32 = 0;
        for result in rdr.deserialize() {
            let record: WeatherRecord = result?;
            jdd = record.date.ordinal() as i32;
            let j = jdd - unsafe { DayStart };
            if j < 0 {
                continue;
            }
            unsafe {
                Clim[j as usize].nDay = jdd;
                // convert $MJ\ m^{-2}$ to langleys
                Clim[j as usize].Rad = record.irradiation * 23.884;
                Clim[j as usize].Tmax = record.tmax;
                Clim[j as usize].Tmin = record.tmin;
                Clim[j as usize].Wind =
                    if self.site.average_wind_speed.is_some() && record.wind.is_none() {
                        self.site.average_wind_speed.unwrap()
                    } else {
                        record.wind.unwrap_or(0.)
                    };
                Clim[j as usize].Tdew = record.tdew.unwrap_or(estimate_dew_point(
                    record.tmax,
                    self.site.estimate_dew_point.0,
                    self.site.estimate_dew_point.1,
                ));
                Clim[j as usize].Rain = record.rain;
            }
        }
        unsafe {
            LastDayWeatherData = jdd;
        }
        let mut idef: usize = 0;
        let mut icult: usize = 0;
        let mut ipix: usize = 0;
        unsafe {
            NumNitApps = 0;
            NumIrrigations = 0;
            NumWaterTableData = 0;
            for i in 0..5 {
                DefoliationDate[i] = 0;
                DefoliationMethod[i] = 0;
                DefoliantAppRate[i] = 0.;
            }

            for ao in &self.agronomy_operations {
                match ao {
                    AgronomyOperation::irrigation {
                        date,
                        amount,
                        predict,
                        method,
                        drip_x,
                        drip_y,
                        max_amount,
                        stop_predict_date,
                    } => {
                        if *predict {
                            MaxIrrigation = max_amount.unwrap();
                            DayStartPredIrrig = date.ordinal() as i32;
                            DayStopPredIrrig = stop_predict_date.unwrap().ordinal() as i32;
                            if let IrrigationMethod::Drip = method {
                                LocationColumnDrip =
                                    utils::slab_horizontal_location(*drip_x, RowSpace)? as i32;
                                LocationLayerDrip = utils::slab_vertical_location(*drip_y)? as i32;
                            }
                            IrrigMethod = *method as i32;
                        } else {
                            Irrig[NumIrrigations as usize].day = date.ordinal() as i32;
                            Irrig[NumIrrigations as usize].amount = *amount;
                            if let IrrigationMethod::Drip = method {
                                Irrig[NumIrrigations as usize].LocationColumnDrip =
                                    utils::slab_horizontal_location(*drip_x, RowSpace)? as i32;
                                Irrig[NumIrrigations as usize].LocationLayerDrip =
                                    utils::slab_vertical_location(*drip_y)? as i32;
                            }
                            Irrig[NumIrrigations as usize].method = *method as i32;
                            NumIrrigations += 1;
                        }
                    }
                    AgronomyOperation::fertilization {
                        date,
                        urea,
                        nitrate,
                        ammonium,
                        method,
                        drip_x,
                        drip_y,
                    } => {
                        NFertilizer[NumNitApps as usize].day = date.ordinal() as i32;
                        NFertilizer[NumNitApps as usize].amtamm = *ammonium;
                        NFertilizer[NumNitApps as usize].amtnit = *nitrate;
                        NFertilizer[NumNitApps as usize].amtura = *urea;
                        match method {
                            FertilizationMethod::Sidedress | FertilizationMethod::Drip => {
                                NFertilizer[NumNitApps as usize].ksdr =
                                    utils::slab_horizontal_location(*drip_x, RowSpace)? as i32;
                                NFertilizer[NumNitApps as usize].lsdr =
                                    utils::slab_vertical_location(*drip_y)? as i32;
                            }
                            _ => {}
                        }
                        NFertilizer[NumNitApps as usize].mthfrt = *method as i32;
                        NumNitApps += 1;
                    }
                    AgronomyOperation::defoliation {
                        date,
                        open_ratio,
                        predict,
                        ppa,
                    } => {
                        DefoliationDate[idef] = date.ordinal() as i32;
                        DefoliantAppRate[idef] = if *predict { -99.9 } else { *ppa };
                        DefoliationMethod[idef] = *open_ratio;
                        DayFirstDef = DefoliationDate[0];
                        idef += 1;
                    }
                    AgronomyOperation::cultivation { date, depth } => {
                        CultivationDate[icult] = date.ordinal() as i32;
                        CultivationDepth[icult] = *depth;
                        icult += 1;
                    }
                    AgronomyOperation::pix { date, method, ppa } => {
                        pixday[ipix] = date.ordinal() as i32;
                        pixmth[ipix] = *method as i32;
                        pixppa[ipix] = *ppa;
                        ipix += 1;
                    }
                    AgronomyOperation::watertable { date, level, ecs } => {
                        DayWaterTableInput[NumWaterTableData as usize] = date.ordinal() as i32;
                        LevelsOfWaterTable[NumWaterTableData as usize] = *level;
                        ElCondSatSoil[NumWaterTableData as usize] = *ecs;
                        NumWaterTableData += 1;
                    }
                }
            }
            if self.light_intercept_method == 2 {
                light_intercept_parameter = 0.;
                for i in 0..20 {
                    light_intercept_parameters[i] = self.light_intercept_parameters.unwrap()[i];
                    light_intercept_parameter += light_intercept_parameters[i];
                }
            }
        }
        self.read_soil_impedance(self.soil_impedance.as_ref().unwrap())?;
        unsafe {
            InitSoil(&self.soil_layers, &self.soil_hydraulic);
            if self.plant_maps.is_some() {
                self.read_plant_map_input();
            }
            InitializeRootData();
            //     initialize some variables at the start of simulation.
            SoilNitrogenAtStart = TotalSoilNo3N + TotalSoilNh4N + TotalSoilUreaN;
            PlantWeightAtStart = TotalRootWeight + TotalStemWeight + TotalLeafWeight() + ReserveC;
        }
        // If this is the first time the function is executed, get the ambient CO2 correction.
        self.ambient_CO2_factor = utils::ambient_CO2_factor(self.start_date.year());
        Ok(())
    }

    /// Write the output CSV file's header.
    pub fn output_file_headers(self: &Self) -> Result<(), Box<dyn std::error::Error>> {
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
    /// Global or file scope variables set:
    /// * [gh2oc]
    /// * [impede]
    /// * [inrim]
    /// * [ncurve]
    /// * [tstbd].
    ///
    fn read_soil_impedance(
        self: &Self,
        path: &std::path::Path,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut rdr = csv::Reader::from_path(path)?;
        unsafe {
            ncurve = 0;
            inrim = 0;
        }
        for result in rdr.records() {
            let record = result?;
            for (i, item) in record.iter().enumerate() {
                if i == 0 {
                    let water_content: f64 = item.parse()?;
                    unsafe {
                        gh2oc[ncurve as usize] = water_content;
                    }
                } else {
                    unsafe {
                        impede[i][ncurve as usize] = item.parse()?;
                    }
                }
            }
            unsafe {
                ncurve += 1;
            }
        }
        for (j, header) in rdr.headers()?.iter().enumerate() {
            if j > 0 {
                unsafe {
                    for i in 0..ncurve as usize {
                        tstbd[j - 1][i] = header.parse()?;
                    }
                }
            }
        }
        Ok(())
    }

    /// This sunbroutine opens and reads an ascii file with input of observed plant map adjustment data. It is used to
    /// adjust the simulation.
    ///
    /// It is called by [Profile::initialize()].
    ///
    /// The following global variables are set:
    /// * [MapDataGreenBollNum]
    /// * [MapDataDate],
    /// * [MapDataMainStemNodes]
    /// * [MapDataPlantHeight]
    /// * [MapDataSquareNum]
    /// * [MapDataAllSiteNum]
    fn read_plant_map_input(&self) {
        for (i, plant_map) in self.plant_maps.as_ref().unwrap().iter().enumerate() {
            unsafe {
                MapDataDate[i] = plant_map.date.ordinal() as i32; // day of year
                MapDataPlantHeight[i] = plant_map.plant_height; // Plant height, cm
                MapDataMainStemNodes[i] = plant_map.main_stem_nodes; // Number of mainstem nodes
                MapDataSquareNum[i] = plant_map.number_of_squares; // Number of squares per plant
                MapDataGreenBollNum[i] = plant_map.number_of_bolls; // Number of green bolls per plant
                MapDataAllSiteNum[i] = plant_map.number_of_nodes; // Number of total sites per plant
            }
        }
    }

    pub fn write_record(self: &Self) -> Result<(), Box<dyn std::error::Error>> {
        let mut f = std::fs::OpenOptions::new()
            .write(true)
            .append(true)
            .open(self.path.parent().unwrap().join("output.csv"))?;
        let mut record = vec![
            unsafe {
                chrono::NaiveDate::from_yo(iyear, Daynum as u32)
                    .format("%F")
                    .to_string()
            },
            unsafe { LightIntercept.to_string() },
            unsafe { LintYield.to_string() },
            unsafe { LeafAreaIndex.to_string() },
            unsafe {
                ((CottonWeightOpenBolls + CottonWeightGreenBolls) * PlantPopulation / 1000.)
                    .to_string()
            },
            unsafe { PlantHeight.to_string() },
            unsafe { NumFruitBranches[0].to_string() },
            unsafe { TotalLeafWeight().to_string() },
            unsafe { TotalPetioleWeight.to_string() },
            unsafe { TotalStemWeight.to_string() },
            unsafe { NumSquares.to_string() },
            unsafe { NumGreenBolls.to_string() },
            unsafe { NumOpenBolls.to_string() },
            unsafe { TotalSquareWeight.to_string() },
            unsafe {
                (CottonWeightOpenBolls
                    + CottonWeightGreenBolls
                    + BurrWeightGreenBolls
                    + BurrWeightOpenBolls)
                    .to_string()
            },
            unsafe { TotalRootWeight.to_string() },
            unsafe {
                (if Daynum >= DayEmerge && isw > 0 {
                    PlantWeight - TotalRootWeight
                } else {
                    0.
                } * PlantPopulation
                    / 1000.)
                    .to_string()
            },
            unsafe { VolWaterContent[3][0].to_string() },
            unsafe { VolWaterContent[5][0].to_string() },
            unsafe { VolWaterContent[7][0].to_string() },
            unsafe { VolWaterContent[3][4].to_string() },
            unsafe { VolWaterContent[5][4].to_string() },
            unsafe { VolWaterContent[7][4].to_string() },
            unsafe { VolWaterContent[3][8].to_string() },
            unsafe { VolWaterContent[5][8].to_string() },
            unsafe { VolWaterContent[7][8].to_string() },
            unsafe { VolWaterContent[3][12].to_string() },
            unsafe { VolWaterContent[5][12].to_string() },
            unsafe { VolWaterContent[7][12].to_string() },
        ];
        unsafe {
            for v in &LeafAreaIndexes {
                record.push(v.to_string());
            }
        }
        writeln!(f, "{}", record.join(","))?;
        Ok(())
    }

    /// This function is called from [Profile::run()].
    /// It checks if plant adjustment data are available for this day and calls the necessary functions to compute
    /// adjustment.
    ///
    /// It calls:
    /// * [PlantAdjustments()]
    /// * [Profile::simulate_this_day()]
    /// * [WriteStateVariables()]
    ///
    /// The following global variable are referenced:
    /// * [DayEmerge]
    /// * [Daynum]
    /// * [Kday]
    ///
    /// The following global variable are set:
    /// * [MapDataDate]
    /// * [nadj]
    /// * [NumAdjustDays]
    pub fn adjust(self: &mut Self) -> Result<bool, Cotton2KError> {
        // Check if plant map data are available for this day. If there are no more adjustments, return.
        let mut sumsad = 0; // sum for checking if any adjustments are due
        for i in 0..30 {
            unsafe {
                sumsad += MapDataDate[i];
            }
        }
        if sumsad <= 0 {
            return Ok(false);
        }
        // Loop for all adjustment data, and check if there is an adjustment for this day.
        for i in 0..30 {
            unsafe {
                if Daynum == MapDataDate[i] {
                    // Compute NumAdjustDays, the number of days for retroactive adjustment. This can not be more than
                    // 12 days, limited by the date of the previous adjustment.
                    NumAdjustDays = Kday - self.kprevadj as i32;
                    if NumAdjustDays > 12 {
                        NumAdjustDays = 12;
                    }
                    // Loop for six possible adjustments. On each iteration call first PlantAdjustments(), which will
                    // assign true to nadj(jj) if adjustment is necessary, and compute the necessary parameters.
                    for jj in 0..5 as usize {
                        adjust::PlantAdjustments(i, jj as i32);
                        //     If adjustment is necessary, rerun the simulation for the
                        //     previous NumAdjustDays (number
                        //  of days) and call WriteStateVariables() to write state
                        //  variables in scratch structure.
                        if nadj[jj] {
                            for _j1 in 0..NumAdjustDays {
                                match self.simulate_this_day() {
                                    Err(e) => {
                                        if e.level == 0 {
                                            println!("{}", e.message);
                                            break;
                                        }
                                    }
                                    _ => {}
                                };
                                if Kday > 0 {
                                    WriteStateVariables(true);
                                }
                            }
                        }
                    }
                    // After finishing this adjustment date, set kprevadj (date of previous adjustment, to be used for
                    // next adjustment), and assign zero to the present msadte, and to array nadj[].
                    self.kprevadj = (MapDataDate[i] - DayEmerge + 1) as u32;
                    MapDataDate[i] = 0;
                    for jj in 0..5 {
                        nadj[jj] = false;
                    }
                }
            }
        }
        return Ok(true);
    }

    /// This function executes all the simulation computations in a day. It is called from [Profile::run()], and
    /// [Profile::adjust()].
    ///
    /// It calls the following functions:
    /// * [Profile::column_shading()]
    /// * [DayClim()]
    /// * [SoilTemperature()]
    /// * [SoilProcedures()]
    /// * [SoilNitrogen()]
    /// * [SoilSum()]
    /// * [PhysiologicalAge()]
    /// * [Profile::pix()]
    /// * [Defoliate()]
    /// * [Profile::stress()]
    /// * [Profile::get_net_photosynthesis()]
    /// * [PlantGrowth()]
    /// * [CottonPhenology()]
    /// * [PlantNitrogen()]
    /// * [CheckDryMatterBal()]
    /// * [PlantNitrogenBal()]
    /// * [SoilNitrogenBal()]
    /// * [SoilNitrogenAverage()]
    ///
    /// The following global variables are referenced here:
    /// * [DayEmerge]
    /// * [DayFinish]
    /// * [DayStart]
    /// * [Kday]
    /// * [LastDayWeatherData]
    /// * [LeafAreaIndex]
    /// * [pixday]
    ///
    /// The following global variables are set here:
    /// * [bEnd]
    /// * [DayInc]
    /// * [Daynum]
    /// * [DayOfSimulation]
    /// * [isw]
    /// * [Kday]
    pub fn simulate_this_day(&mut self) -> Result<(), Cotton2KError> {
        unsafe {
            // Compute Daynum (day of year), Date, and DayOfSimulation (days from start of simulation).
            Daynum += 1;
            DayOfSimulation = Daynum - DayStart + 1;
            //    Compute Kday (days from emergence).
            if DayEmerge <= 0 {
                Kday = 0;
            } else {
                Kday = Daynum - DayEmerge + 1;
            }
            if Kday < 0 {
                Kday = 0;
            }
            // The following functions are executed each day (also before emergence).
            self.column_shading(); // computes light interception and soil shading.
            DayClim(); // computes climate variables for today.
            SoilTemperature(); // executes all modules of soil and canopy temperature.
            SoilProcedures()?; // executes all other soil processes.
            SoilNitrogen(); // computes nitrogen transformations in the soil.
            SoilSum(); // computes totals of water and N in the soil.

            // The following is executed each day after plant emergence:
            if Daynum >= DayEmerge && isw > 0 {
                // If this day is after emergence, assign to isw the value of 2.
                isw = 2;
                DayInc = PhysiologicalAge(); // computes physiological age
                if pixday[0] > 0 {
                    self.pix(); // effects of pix applied.
                }
                Defoliate(); // effects of defoliants applied.
                self.stress(); // computes water stress factors.
                self.get_net_photosynthesis()?; // computes net photosynthesis.
                PlantGrowth(); // executes all modules of plant growth.
                CottonPhenology(); // executes all modules of plant phenology.
                PlantNitrogen(); // computes plant nitrogen allocation.
                CheckDryMatterBal(); // checks plant dry matter balance.

                // If the relevant output flag is not zero, compute soil nitrogen balance and soil nitrogen averages by
                // layer, and write this information to files.
                if false {
                    PlantNitrogenBal(); // checks plant nitrogen balance.
                    SoilNitrogenBal(); // checks soil nitrogen balance.
                    SoilNitrogenAverage(); // computes average soil nitrogen by layers.
                }
            }
            // Check if the date to stop simulation has been reached, or if this is the last day with available weather
            // data. Simulation will also stop when no leaves remain on the plant.
            if Daynum >= LastDayWeatherData {
                return Err(Cotton2KError {
                    level: 0,
                    message: String::from("No more weather data!"),
                });
            }
            if Kday > 10 && LeafAreaIndex < 0.0002 {
                return Err(Cotton2KError {
                    level: 0,
                    message: String::from("Leaf area index is too small!"),
                });
            }
        }
        Ok(())
    }

    /// This function computes light interception by crop canopy and shading of soil columns by the plants. It is called from SimulateThisDay().
    ///
    /// The following global variables are referenced here:
    /// * [DayEmerge]
    /// * [Daynum]
    /// * [isw]
    /// * [LeafAreaIndex]
    /// * [PlantHeight]
    /// * [PlantRowColumn]
    /// * [nk]
    /// * [RowSpace]
    ///
    /// The following global variables are set here:
    /// * [LightIntercept]
    /// * [rracol]
    fn column_shading(&mut self) {
        unsafe {
            // Before emergence: no light interception and no shading.
            // LightIntercept is assigned zero, and the rracol array is assigned 1.
            if Daynum < DayEmerge || isw <= 0 || DayEmerge <= 0 {
                LightIntercept = 0.;
                for k in 0..nk as usize {
                    rracol[k] = 1.;
                }
                return;
            }
            // Compute the maximum leaf area index until this day (lmax).
            if Daynum <= DayEmerge {
                self.lmax = 0.;
            } else if LeafAreaIndex > self.lmax {
                self.lmax = LeafAreaIndex;
            }
            // Light interception is computed by two methods:
            //
            // 1. It is assumed to be proportional to the ratio of plant height to
            //    row spacing.
        }
        // light interception computed from plant height.
        let zint = unsafe { 1.0756 * PlantHeight / RowSpace };
        unsafe {
            if light_intercept_method == LIGHT_INTERCEPT_METHOD_LIGHT_INTERCEPT_METHOD_LAYERED {
                for i in 0..20 {
                    LeafArea[i] = 0.;
                    AverageLeafAge[i] = 0.;
                }
                for i in 0..9 {
                    LeafArea[NodeLayerPreFru[i] as usize] += LeafAreaPreFru[i];
                    AverageLeafAge[NodeLayerPreFru[i] as usize] +=
                        LeafAreaPreFru[i] * AgeOfPreFruNode[i];
                }
                for k in 0..NumVegBranches as usize {
                    for l in 0..NumFruitBranches[k] as usize {
                        LeafArea[NodeLayer[k][l] as usize] += LeafAreaMainStem[k][l];
                        AverageLeafAge[NodeLayer[k][l] as usize] +=
                            LeafAreaMainStem[k][l] * LeafAge[k][l][0];
                        for m in 0..NumNodes[k][l] as usize {
                            LeafArea[NodeLayer[k][l] as usize] += LeafAreaNodes[k][l][m];
                            AverageLeafAge[NodeLayer[k][l] as usize] +=
                                LeafAreaNodes[k][l][m] * LeafAge[k][l][m];
                        }
                    }
                }
                if FirstSquare <= 0 {
                    LeafArea[0] += 0.20 * 0.6;
                }
                let mut light_through = 0.;
                for i in 0..20 {
                    AverageLeafAge[i] /= LeafArea[i];
                    LeafAreaIndexes[i] = LeafArea[i] / PerPlantArea;
                    LightInterceptLayer[i] =
                        1. - (light_intercept_parameters[i] * LeafAreaIndexes[i]).exp();
                    light_through += light_intercept_parameters[i] * LeafAreaIndexes[i];
                }
                LightIntercept = 1. - light_through.exp();
            } else if light_intercept_method
                == LIGHT_INTERCEPT_METHOD_LIGHT_INTERCEPT_METHOD_FRY1980
            {
                LightIntercept = 0.39 * LeafAreaIndex.powf(0.68);
            } else {
                // 2. It is computed as a function of leaf area index. If LeafAreaIndex is not greater than 0.5 lfint is a
                //    linear function of it.

                // light interception computed from leaf area index.
                let lfint = if LeafAreaIndex <= 0.5 {
                    0.80 * LeafAreaIndex
                } else {
                    // If the leaf area index is greater than 0.5, lfint is computed as an exponential function of
                    // LeafAreaIndex.
                    1. - (0.07 - 1.16 * LeafAreaIndex).exp()
                };
                // If lfint is greater then zint, LightIntercept is their average value.
                // Otherwise, if the LeafAreaIndex is decreasing, it is lfint. Else it is zint.
                LightIntercept = if lfint > zint {
                    0.5 * (zint + lfint)
                } else if LeafAreaIndex < self.lmax {
                    lfint
                } else {
                    zint
                };
            }
            // The value of LightIntercept is between zero and one.
            if LightIntercept < 0. {
                LightIntercept = 0.;
            }
            if LightIntercept > 1. {
                LightIntercept = 1.;
            }
        }
        // Loop of soil columns.
        let mut sw = 0.; // sum of column widths
        let mut sw0: f64 = 0.; // sum of column widths up to location of plant row.
        let mut sw1: f64; // distance from middle of a column to the plant row, cm.
        let mut j;
        let mut k0; // number of columns from plant row location.
        unsafe {
            for k in 0..nk as usize {
                if k <= PlantRowColumn as usize {
                    // When the column is on the left of the plant row.
                    j = (PlantRowColumn - k as i32) as usize;
                    sw += wk[j];
                    sw0 = sw;
                    sw1 = sw - wk[j] / 2.;
                    k0 = j;
                } else {
                    // When the column is on the right of the plant row.
                    sw += wk[k];
                    sw1 = sw - sw0 - wk[k] / 2.;
                    k0 = k;
                }

                // Relative shading is computed as a function of sw1 and plant height, modified by LightIntercept.
                // [rracol] is the fraction of radiation received by a soil column. Iys minimum value is 0.05 .

                // relative shading of a soil column.
                let shade = if sw1 >= PlantHeight {
                    0.
                } else {
                    let mut result = 1. - (sw1 / PlantHeight).powi(2);
                    if LightIntercept < zint && LeafAreaIndex < self.lmax {
                        result *= LightIntercept / zint;
                    }
                    result
                };
                rracol[k0] = 1. - shade;
                if rracol[k0] < 0.05 {
                    rracol[k0] = 0.05;
                }
            }
        }
    }

    /// effects of pix applied.
    ///
    /// TODO
    fn pix(self: &Self) {}

    /// This function simulates the net photosynthesis of cotton  plants. It is called daily by
    /// [Profile::simulate_this_day()]. This is essentially the routine of GOSSYM with minor changes.
    ///
    /// The following global and file scope variables are referenced here:
    /// * [BurrWeightOpenBolls]
    /// * [CottonWeightOpenBolls]
    /// * [DayLength]
    /// * [Daynum]
    /// * [DayEmerge]
    /// * [DayTimeTemp]
    /// * [iyear]
    /// * [Kday]
    /// * [LeafNConc]
    /// * [LightIntercept]
    /// * [PerPlantArea]
    /// * [PlantWeight]
    /// * [Profile::ptsred]
    /// * [StemWeight]
    /// * [TotalLeafWeight()]
    ///
    /// The following global variables are set here:
    /// * [bEnd]
    /// * [CumNetPhotosynth]
    /// * [NetPhotosynthesis]
    ///
    /// References:
    ///
    /// * Baker et. al. (1972). Simulation of Growth and Yield in Cotton: I. Gross photosynthesis, respiration and
    ///   growth. Crop Sci. 12:431-435.
    /// * Harper et. al. (1973) Carbon dioxide and the photosynthesis of field crops. A metered carbon dioxide release
    ///   in cotton under field conditions.  Agron. J. 65:7-11.
    /// * Baker (1965) Effects of certain environmental factors on net assimilation in cotton. Crop Sci. 5:53-56 (Fig 5).
    unsafe fn get_net_photosynthesis(&mut self) -> Result<(), Cotton2KError> {
        //  constants:
        const gsubr: f64 = 0.375; // the growth resiration factor.
        const rsubo: f64 = 0.0032; // maintenance respiration factor.
        const vpnet: [f64; 4] = [1.30, 0.034, 0.010, 0.32];
        const co2parm: [f64; 45] =
            // parameters used to correct photosynthesis for ambient CO2 concentration.
            [
                1.0235, 1.0264, 1.0285, 1.0321, 1.0335, 1.0353, 1.0385, 1.0403, 1.0431, 1.0485,
                1.0538, 1.0595, 1.0627, 1.0663, 1.0716, 1.0752, 1.0784, 1.0823, 1.0880, 1.0923,
                1.0968, 1.1019, 1.1087, 1.1172, 1.1208, 1.1243, 1.1311, 1.1379, 1.1435, 1.1490,
                1.1545, 1.1601, 1.1656, 1.1712, 1.1767, 1.1823, 1.1878, 1.1934, 1.1990, 1.2045,
                1.2101, 1.2156, 1.2212, 1.2267, 1.2323,
            ];
        // Note: co2parm is for icrease in ambient CO2 concentration changes from 1959 (308 ppm).
        // The first 28 values (up to 1987) are from GOSSYM. The other values (up to 2004) are derived from data of the
        // Carbon Dioxide Information Analysis Center (CDIAC).
        //
        // Exit the function and end simulation if there are no leaves.

        if TotalLeafWeight() <= 0. {
            return Err(Cotton2KError {
                level: 0,
                message: String::from("Leaf weight is less than 0!"),
            });
        }

        // Get the CO2 correction factor (pnetcor) for photosynthesis, using AmbientCO2Factor and a factor that may be
        // variety specific (vpnet[0]).

        // correction factor for gross photosynthesis.
        let pnetcor =
            self.ambient_CO2_factor * vpnet[0] * self.co2_enrichment.unwrap_or_default().factor;
        // Compute ptnfac, the effect of leaf N concentration on photosynthesis, using an empirical relationship.
        // correction factor for low nitrogen content in leaves.
        let mut ptnfac =
            vpnet[3] + (LeafNConc - vpnet[2]) * (1. - vpnet[3]) / (vpnet[1] - vpnet[2]);
        if ptnfac > 1. {
            ptnfac = 1.;
        }
        if ptnfac < vpnet[3] {
            ptnfac = vpnet[3];
        }
        // Convert the average daily short wave radiation from langley per day, to Watts per square meter (wattsm).
        // average daily global radiation, W m-2.
        let wattsm =
            GetFromClim(CLIMATE_METRIC_CLIMATE_METRIC_IRRD, Daynum) * 697.45 / (DayLength * 60.);
        // Compute pstand as an empirical function of wattsm (based on Baker et al., 1972).
        // gross photosynthesis for a non-stressed full canopy.
        let pstand = 2.3908 + wattsm * (1.37379 - wattsm * 0.00054136);
        // Convert it to gross photosynthesis per plant (pplant), using PerPlantArea and corrections for light
        // interception by canopy, ambient CO2 concentration, water stress and low N in the leaves.
        let mut pplant = 0.;
        // actual gross photosynthetic rate, g per plant per day.
        if light_intercept_method == LIGHT_INTERCEPT_METHOD_LIGHT_INTERCEPT_METHOD_LAYERED {
            let mut pstand_remain = pstand;
            for i in (0..20).rev() {
                if pstand_remain <= 0. {
                    break;
                }
                if LightInterceptLayer[i] <= 0. {
                    continue;
                }
                let page = 1. - (AverageLeafAge[i] / drop_leaf_age(LeafArea[i])).powi(2);
                let mut pplant_inc = 0.001
                    * pstand_remain
                    * LightInterceptLayer[i]
                    * PerPlantArea
                    * self.ptsred
                    * pnetcor
                    * ptnfac
                    * page;
                if pplant_inc > pstand_remain {
                    pplant_inc = pstand_remain;
                }
                pplant += pplant_inc;
                pstand_remain -= pplant_inc;
            }
        } else {
            pplant =
                0.001 * pstand * LightIntercept * PerPlantArea * self.ptsred * pnetcor * ptnfac;
        };
        // Compute the photorespiration factor (rsubl) as a linear function af average day time temperature.
        let rsubl = 0.0032125 + 0.0066875 * DayTimeTemp; // photorespiration factor.

        // Photorespiration (lytres) is computed as a proportion of gross photosynthetic rate.
        let lytres = rsubl * pplant; // rate of photorespiration, g per plant per day.

        // Old stems are those more than voldstm = 32 calendar days old.
        // Maintenance respiration is computed on the basis of plant dry weight, minus the old stems and the dry tissue
        // of opened bolls.
        let voldstm = 32;
        let kkday = Kday - voldstm; // day of least recent actively growing stems.

        // weight of old stems.
        let oldstmwt = if kkday < 1 {
            0.
        } else {
            StemWeight[kkday as usize]
        };
        // maintenance respiration, g per plant per day.
        let bmain = (PlantWeight - CottonWeightOpenBolls - BurrWeightOpenBolls - oldstmwt) * rsubo;
        // Net photosynthesis is computed by substracting photo-respiration and maintenance respiration from the gross
        // rate of photosynthesis. To avoid computational problems, make sure that pts is positive and non-zero.
        // intermediate computation of NetPhotosynthesis.
        let mut pts = pplant - lytres - bmain;
        if pts < 0.00001 {
            pts = 0.00001;
        }
        // The growth respiration (gsubr) supplies energy for converting the supplied carbohydrates to plant tissue dry
        // matter. 0.68182 converts CO2 to CH2O. NetPhotosynthesis is the computed net photosynthesis, in g per plant
        // per day.

        NetPhotosynthesis = pts / (1. + gsubr) * 0.68182;
        // CumNetPhotosynth is the cumulative value of NetPhotosynthesis, from day of emergence.
        CumNetPhotosynth += NetPhotosynthesis;
        Ok(())
    }

    /// This function computes the water stress variables affecting the cotton plants.
    /// It is called by [Profile::simulate_this_day()] and calls [LeafWaterPotential()].
    ///
    /// The following global variables are referenced here:
    /// * [Kday]
    /// * [LwpMin]
    /// * [LwpMax]
    /// The following global variables are set here:
    /// * [AverageLwp]
    /// * [AverageLwpMin]
    /// * [LwpMinX]
    /// * [LwpX]
    /// * [Profile::ptsred]
    /// * [WaterStress]
    /// * [WaterStressStem]
    unsafe fn stress(&mut self) {
        // The following constant parameters are used:
        const vstrs: [f64; 9] = [-3.0, 3.229, 1.907, 0.321, -0.10, 1.230, 0.340, 0.30, 0.05];
        // Call LeafWaterPotential() to compute leaf water potentials.
        LeafWaterPotential();
        // The running averages, for the last three days, are computed: AverageLwpMin is the average of LwpMin, and
        // AverageLwp of LwpMin + LwpMax.
        AverageLwpMin += (LwpMin - LwpMinX[2]) / 3.;
        AverageLwp += (LwpMin + LwpMax - LwpX[2]) / 3.;
        for i in [2, 1] {
            LwpMinX[i] = LwpMinX[i - 1];
            LwpX[i] = LwpX[i - 1];
        }
        LwpMinX[0] = LwpMin;
        LwpX[0] = LwpMin + LwpMax;
        //     No stress effects before 5th day after emergence.
        if Kday < 5 {
            self.ptsred = 1.;
            WaterStress = 1.;
            WaterStressStem = 1.;
            return;
        }
        // The computation of ptsred, the effect of moisture stress on the photosynthetic rate, is based on the
        // following work:
        // * Ephrath, J.E., Marani, A., Bravdo, B.A., 1990. Effects of moisture stress on stomatal resistance and
        //   photosynthetic rate in cotton (Gossypium hirsutum) 1. Controlled levels of stress. Field Crops Res.
        //   23:117-131.
        //
        // It is a function of AverageLwpMin (average LwpMin for the last three days).
        if AverageLwpMin < vstrs[0] {
            AverageLwpMin = vstrs[0];
        }
        self.ptsred = vstrs[1] + AverageLwpMin * (vstrs[2] + vstrs[3] * AverageLwpMin);
        if self.ptsred > 1. {
            self.ptsred = 1.;
        }
        // The general moisture stress factor (WaterStress) is computed as an empirical function of AverageLwp. psilim,
        // the value of AverageLwp at the maximum value of the function, is used for truncating it.

        // The minimum value of WaterStress is 0.05, and the maximum is 1.
        let psilim = -0.5 * vstrs[5] / vstrs[6]; // limiting value of AverageLwp.
        if AverageLwp > psilim {
            WaterStress = 1.;
        } else {
            WaterStress = vstrs[4] - AverageLwp * (vstrs[5] + vstrs[6] * AverageLwp);
            if WaterStress > 1. {
                WaterStress = 1.;
            }
            if WaterStress < 0.05 {
                WaterStress = 0.05;
            }
        }
        // Water stress affecting plant height and stem growth (WaterStressStem) is assumed to be more severe than
        // WaterStress, especially at low WaterStress values.
        WaterStressStem = WaterStress * (1. + vstrs[7] * (2. - WaterStress)) - vstrs[7];
        if WaterStressStem < vstrs[8] {
            WaterStressStem = vstrs[8];
        }
    }
}

fn fmin(a: f64, b: f64) -> f64 {
    if a < b {
        a
    } else {
        b
    }
}
fn fmax(a: f64, b: f64) -> f64 {
    if a > b {
        a
    } else {
        b
    }
}

/// This function simulates the leaf water potential of cotton plants. It has been adapted from the model of Moshe Meron (The relation of cotton leaf water potential to soil water content in the irrigated management range. PhD dissertation, UC Davis, 1984).
///
///  It is called from [Profile::stress()]. It calls [wcond()] and [LeafResistance()].
///
/// The following global variables are referenced here:
/// * [AgeOfPreFruNode]
/// * [AverageSoilPsi]
/// * [beta]
/// * [dl]
/// * [Kday]
/// * [LeafAge]
/// * [NumFruitBranches]
/// * [NumLayersWithRoots]
/// * [NumNodes]
/// * [NumPreFruNodes]
/// * [NumVegBranches]
/// * [PlantHeight]
/// * [PoreSpace]
/// * [ReferenceETP]
/// * [RootColNumLeft]
/// * [RootColNumRight]
/// * [RootWtCapblUptake]
/// * [SaturatedHydCond]
/// * [SoilPsi]
/// * [thad]
/// * [thts]
/// * [VolWaterContent]
/// * [wk]
///
/// The following global variables are set here:
/// * [LwpMin]
/// * [LwpMax]
fn LeafWaterPotential() {
    // Constant parameters used:
    const cmg: f64 = 3200.; // length in cm per g dry weight of roots, based on an average root diameter of 0.06 cm, and a specific weight of 0.11 g  dw per cubic cm.
    const psild0: f64 = -1.32; // maximum values of LwpMin
    const psiln0: f64 = -0.40; // maximum values of LwpMax.
    const rtdiam: f64 = 0.06; // average root diameter in cm.
    const vpsil: [f64; 13] = [
        0.48, -5.0, 27000., 4000., 9200., 920., 0.000012, -0.15, -1.70, -3.5, 0.1e-9, 0.025, 0.80,
    ];
    // Leaf water potential is not computed during 10 days after emergence. Constant values are assumed for this period.
    unsafe {
        if Kday <= 10 {
            LwpMax = psiln0;
            LwpMin = psild0;
            return;
        }
    }
    // Compute shoot resistance (rshoot) as a function of plant height.
    let rshoot: f64 = vpsil[0] * unsafe { PlantHeight } / 100.; // shoot resistance, Mpa hours per cm.

    // Assign zero to summation variables
    let mut psinum = 0f64; // sum of RootWtCapblUptake for all soil cells with roots.
    let mut rootvol = 0f64; // sum of volume of all soil cells with roots.
    let mut rrlsum = 0f64; // weighted sum of reciprocals of rrl.
    let rroot; // root resistance, Mpa hours per cm.
    let mut sumlv = 0f64; // weighted sum of root length, cm, for all soil cells with roots.
    let mut vh2sum = 0f64; // weighted sum of soil water content, for all soil cells with roots.

    // Loop over all soil cells with roots. Check if RootWtCapblUptake is greater than vpsil[10].
    // All average values computed for the root zone, are weighted by RootWtCapblUptake (root weight capable of uptake), but the weight assigned will not be greater than vpsil[11].
    unsafe {
        for l in 0..NumLayersWithRoots as usize {
            for k in (RootColNumLeft[l] as usize)..RootColNumLeft[l] as usize {
                if RootWtCapblUptake[l][k] >= vpsil[10] {
                    psinum += fmin(RootWtCapblUptake[l][k], vpsil[11]);
                    sumlv += fmin(RootWtCapblUptake[l][k], vpsil[11]) * cmg;
                    rootvol += dl[l] * wk[k];
                    // root resistance per g of active roots.
                    let rrl = if SoilPsi[l][k] <= vpsil[1] {
                        vpsil[2] / cmg
                    } else {
                        (vpsil[3] - SoilPsi[l][k] * (vpsil[4] + vpsil[5] * SoilPsi[l][k])) / cmg
                    };
                    rrlsum += fmin(RootWtCapblUptake[l][k], vpsil[11]) / rrl;
                    vh2sum += VolWaterContent[l][k] * fmin(RootWtCapblUptake[l][k], vpsil[11]);
                }
            }
        }
    }
    // Compute average root resistance (rroot) and average soil water content (vh2).
    let dumyrs: f64; // intermediate variable for computing cond.
    let vh2: f64; // average of soil water content, for all soil soil cells with roots.
    if psinum > 0. && sumlv > 0. {
        rroot = psinum / rrlsum;
        vh2 = vh2sum / psinum;
        dumyrs = fmax(
            1.001,
            (1. / (std::f64::consts::PI * sumlv / rootvol)).sqrt() / rtdiam,
        );
    } else {
        rroot = 0.;
        vh2 = unsafe { thad[0] };
        dumyrs = 1.001;
    }
    // Compute hydraulic conductivity (cond), and soil resistance near the root surface (rsoil). soil hydraulic conductivity near the root surface.
    let cond = fmax(
        unsafe {
            wcond(
                vh2,
                thad[0],
                thts[0],
                beta[0],
                SaturatedHydCond[0],
                PoreSpace[0],
            )
        } / 24.
            * 2.
            * sumlv
            / rootvol
            / dumyrs.ln(),
        vpsil[6],
    );
    let rsoil = 0.0001 / (2. * std::f64::consts::PI * cond); // soil resistance, Mpa hours per cm.

    // Compute leaf resistance (LeafResistance) as the average of the resistances of all existing leaves. The resistance of an individual leaf is a function of its age. Function LeafResistance is called to compute it.
    // This is executed for all the leaves of the plant.
    let mut numl = 0; // number of leaves.
    let mut sumrl = 0f64; // sum of leaf resistances for all the plant.
    unsafe {
        for j in 0..NumPreFruNodes as usize
        // loop prefruiting nodes
        {
            numl += 1;
            sumrl += LeafResistance(AgeOfPreFruNode[j]);
        }
    }
    //
    let mut nbrch: i32; // number of fruiting branches on a vegetative branch.
    let mut nnid: i32; // number of nodes on a fruiting branch.
    unsafe {
        for k in 0..NumVegBranches as usize {
            // loop for all other nodes
            nbrch = NumFruitBranches[k];
            for l in 0..nbrch as usize {
                nnid = NumNodes[k][l];
                for m in 0..nnid as usize {
                    numl += 1;
                    sumrl += LeafResistance(LeafAge[k][l][m]);
                }
            }
        }
    }
    let rleaf = sumrl / numl as f64; // leaf resistance, Mpa hours per cm.

    // The total resistance to transpiration, MPa hours per cm, (rtotal) is computed.
    let rtotal = rsoil + rroot + rshoot + rleaf;
    // Compute maximum (early morning) leaf water potential, LwpMax, from soil water potential (AverageSoilPsi, converted from bars to MPa).
    // Check for minimum and maximum values.
    unsafe {
        LwpMax = vpsil[7] + 0.1 * AverageSoilPsi;
        if LwpMax < vpsil[8] {
            LwpMax = vpsil[8];
        }
        if LwpMax > psiln0 {
            LwpMax = psiln0;
        }
    }
    //     Compute minimum (at time of maximum transpiration rate) leaf water
    //     potential, LwpMin, from
    //  maximum transpiration rate (etmax) and total resistance to transpiration
    //  (rtotal).
    let mut etmax = 0f64; // the maximum hourly rate of evapotranspiration for this day.
    for ihr in 0..24 {
        //  hourly loop
        unsafe {
            if ReferenceETP[ihr] > etmax {
                etmax = ReferenceETP[ihr];
            }
        }
    }
    unsafe {
        LwpMin = LwpMax - 0.1 * fmax(etmax, vpsil[12]) * rtotal;
        //     Check for minimum and maximum values.
        if LwpMin < vpsil[9] {
            LwpMin = vpsil[9];
        }
        if LwpMin > psild0 {
            LwpMin = psild0;
        }
    }
}

/// This function manages all the soil related processes, and is executed once each day.
///
/// It is called from [Profile::simulate_this_day()] and it calls the following functions:
/// * [ApplyFertilizer()]
/// * [AveragePsi()]
/// * [CapillaryFlow()]
/// * [ComputeIrrigation()]
/// * [DripFlow()]
/// * [GravityFlow()]
/// * [RootsCapableOfUptake()]
/// * [WaterUptake()]
/// * [WaterTable()]
///
/// The following global variables are referenced here:
/// * [ActualTranspiration]
/// * [Clim]
/// * [DayEmerge]
/// * [Daynum]
/// * [DayStartPredIrrig]
/// * [DayStopPredIrrig]
/// * [dl]
/// * [Irrig]
/// * [IrrigMethod]
/// * [isw]
/// * [Kday]
/// * [MaxIrrigation]
/// * [nk]
/// * [nl]
/// * [NumIrrigations]
/// * [NumWaterTableData]
/// * [PerPlantArea]
/// * [SoilPsi]
/// * [RowSpace]
/// * [SupplyNH4N]
/// * [SupplyNO3N]
/// * [VolWaterContent]
/// * [wk]
///
/// The following global variables are set here:
/// * [AverageSoilPsi]
/// * [CumNitrogenUptake]
/// * [CumTranspiration]
/// * [CumWaterAdded]
/// * [LocationColumnDrip]
/// * [LocationLayerDrip]
/// * [noitr]
fn SoilProcedures() -> Result<(), Cotton2KError> {
    // The following constant parameters are used:
    const cpardrip: f64 = 0.2;
    const cparelse: f64 = 0.4;
    // Call function ApplyFertilizer() for nitrogen fertilizer application.
    unsafe {
        ApplyFertilizer();
    }
    let mut DripWaterAmount = 0f64; // amount of water applied by drip irrigation

    // amount of water applied by non-drip irrigation or rainfall
    // Check if there is rain on this day
    let mut WaterToApply = unsafe { GetFromClim(CLIMATE_METRIC_CLIMATE_METRIC_RAIN, Daynum) };
    // If irrigation is to be predicted for this day, call ComputeIrrigation() to compute the actual amount of irrigation.
    unsafe {
        if MaxIrrigation > 0. {
            if Daynum >= DayStartPredIrrig && Daynum < DayStopPredIrrig {
                ComputeIrrigation();
                if IrrigMethod == 2 {
                    DripWaterAmount = AppliedWater;
                } else {
                    WaterToApply += AppliedWater;
                }
                AppliedWater = 0.;
            }
        }
    }
    // When water is added by an irrigation defined in the input: update the amount of applied water.
    unsafe {
        for i in 0..NumIrrigations as usize {
            if Daynum == Irrig[i].day {
                if Irrig[i].method == 2 {
                    DripWaterAmount += Irrig[i].amount;
                    LocationColumnDrip = Irrig[i].LocationColumnDrip;
                    LocationLayerDrip = Irrig[i].LocationLayerDrip;
                } else {
                    WaterToApply += Irrig[i].amount;
                }
                break;
            }
        }
    }
    unsafe {
        CumWaterAdded += WaterToApply + DripWaterAmount;
    }
    unsafe {
        if Kday > 0 {
            Scratch21[(DayOfSimulation - 1) as usize].amitri = WaterToApply + DripWaterAmount;
        }
    }
    unsafe {
        // The following will be executed only after plant emergence
        if Daynum >= DayEmerge && isw > 0 {
            RootsCapableOfUptake(); // function computes roots capable of uptake for each soil cell
            AverageSoilPsi = AveragePsi(); // function computes the average matric soil water potential in the root zone, weighted by the roots-capable-of-uptake.
            WaterUptake(); // function  computes water and nitrogen uptake by plants.

            // Update the cumulative sums of actual transpiration (CumTranspiration, mm) and total uptake of nitrogen (CumNitrogenUptake, mg N per slab, converted from total N supply, g per plant).
            CumTranspiration += ActualTranspiration;
            CumNitrogenUptake += (SupplyNO3N + SupplyNH4N) * 10. * RowSpace / PerPlantArea;
        }
    }
    // Call function WaterTable() for saturating soil below water table.
    unsafe {
        if NumWaterTableData > 0 {
            WaterTable();
        }
    }
    if WaterToApply > 0. {
        // For rain or surface irrigation.
        // The number of iterations is computed from the thickness of the first soil layer.
        unsafe {
            noitr = (cparelse * WaterToApply / (dl[0] + 2.) + 1.) as i32;
        }
        // the amount of water applied, mm per iteration.
        let applywat = WaterToApply / unsafe { noitr } as f64;
        // The following subroutines are called noitr times per day:
        // If water is applied, GravityFlow() is called when the method of irrigation is not by drippers, followed by CapillaryFlow().
        for _ in 0..unsafe { noitr } as usize {
            unsafe {
                GravityFlow(applywat);
                CapillaryFlow();
            }
        }
    }
    if DripWaterAmount > 0. {
        // For drip irrigation.
        // The number of iterations is computed from the volume of the soil cell in which the water is applied.
        unsafe {
            noitr = (cpardrip * DripWaterAmount
                / (dl[LocationLayerDrip as usize] * wk[LocationColumnDrip as usize])
                + 1.) as i32;
        }
        // the amount of water applied, mm per iteration.
        let applywat = DripWaterAmount / unsafe { noitr } as f64;
        // If water is applied, DripFlow() is called followed by
        // CapillaryFlow().
        for _ in 0..unsafe { noitr } as usize {
            unsafe {
                DripFlow(applywat)?;
                CapillaryFlow();
            }
        }
    }
    // When no water is added, there is only one iteration in this day.
    if WaterToApply + DripWaterAmount <= 0. {
        unsafe {
            noitr = 1;
            CapillaryFlow();
        }
    }
    Ok(())
}
//     This function sets the water saturation of the soil layers below the
//     water
//  table, if it has been defined in the input. It is called from
//  SoilProcedures() if water table data have been input.
//
//     The following global variables are referenced here:
//       Daynum. dl, DayWaterTableInput, ElCondSatSoil, LevelsOfWaterTable,
//       MaxWaterCapacity, nk, nl, NumWaterTableData, PoreSpace,
//       MaxWaterCapacity, RowSpace, wk.
//
//     The following global variables are set here:
//       addwtbl, ElCondSatSoilToday, WaterTableLayer, VolWaterContent.
//
fn WaterTable() {
    if unsafe { NumWaterTableData } <= 0 {
        return;
    }
    //     Find the depth of water table for this day.
    let mut lwtable = 201f64; // level of water table on this day, cm
    unsafe {
        for i in 0..NumWaterTableData as usize {
            if Daynum >= DayWaterTableInput[i] {
                lwtable = LevelsOfWaterTable[i];
                ElCondSatSoilToday = ElCondSatSoil[i];
            }
        }
    }
    //     Find the number of the uppermost layer of water table
    if lwtable > 200. {
        unsafe {
            WaterTableLayer = 1000;
        }
    } else {
        let mut sumdl = 0f64; // sum of depth of consecutive soil layers
        for l in 0..unsafe { nl } as usize {
            sumdl += unsafe { dl[l] };
            if sumdl > lwtable {
                unsafe {
                    WaterTableLayer = l as i32;
                }
                break;
            }
        }
    }
    // The total water entering the soil slab (addwtbl) is computed. It is used to check the water balance in the soil.

    unsafe {
        for l in 0..nl as usize {
            if l as i32 >= WaterTableLayer {
                for k in 0..nk as usize {
                    // previous water content of a cell
                    let vh2ocx = VolWaterContent[l][k];
                    VolWaterContent[l][k] = PoreSpace[l];
                    addwtbl += 10. * (VolWaterContent[l][k] - vh2ocx) * dl[l] * wk[k] / RowSpace;
                }
            } else {
                // Make sure that (in case water table was lowered) water content is not higher than MaxWaterCapacity and adjust addwtbl.
                for k in 0..nk as usize {
                    if VolWaterContent[l][k] > MaxWaterCapacity[l] {
                        // previous water content of a cell
                        let vh2ocx = VolWaterContent[l][k];
                        VolWaterContent[l][k] = MaxWaterCapacity[l];
                        addwtbl +=
                            10. * (VolWaterContent[l][k] - vh2ocx) * dl[l] * wk[k] / RowSpace;
                    }
                }
            }
        }
    }
}

/// This function computes the water redistribution in the soil after irrigation by a drip system.
/// It also computes the resulting redistribution of nitrate and urea N.
/// It is called by SoilProcedures() noitr times per day.
/// It calls function CellDistrib().
///
/// The following argument is used:
/// Drip - amount of irrigation applied by the drip method, mm.
///
/// The following global variables are referenced:
/// * [dl]
/// * [LocationColumnDrip]
/// * [LocationLayerDrip]
/// * [MaxWaterCapacity]
/// * [nk]
/// * [nl]
/// * [NO3FlowFraction]
/// * [PoreSpace]
/// * [RowSpace]
/// * [wk]
///
/// The following global variables are set:
/// * [CumWaterDrained]
/// * [SoilNitrogenLoss]
/// * [VolWaterContent]
/// * [VolNo3NContent]
/// * [VolUreaNContent]
fn DripFlow(Drip: f64) -> Result<(), Cotton2KError> {
    let mut dripw: [f64; 40] = [0f64; 40]; // amount of water applied, or going from one ring of
                                           // soil cells to the next one, cm3. (array)
    let mut dripn: [f64; 40] = [0f64; 40]; // amount of nitrate N applied, or going from one ring of
                                           // soil soil cells to the next one, mg. (array)
    let mut dripu: [f64; 40] = [0f64; 40]; // amount of urea N applied, or going from one ring of
                                           // soil soil cells to the next one, mg. (array)

    // Incoming flow of water (Drip, in mm) is converted to dripw(0), in cm^3 per slab.
    dripw[0] = Drip * unsafe { RowSpace } * 0.10;
    // Wetting the cell in which the emitter is located.
    let mut h2odef; // the difference between the maximum water capacity (at a water content of uplimit) of this ring of soil cell, and the actual water content, cm3.
    let l0 = unsafe { LocationLayerDrip } as usize; //  layer where the drip emitter is situated
    let k0 = unsafe { LocationColumnDrip } as usize; //  column where the drip emitter is situated

    // It is assumed that wetting cannot exceed MaxWaterCapacity of this cell.
    // Compute h2odef, the amount of water needed to saturate this cell.
    h2odef = unsafe { (MaxWaterCapacity[l0] - VolWaterContent[l0][k0]) * dl[l0] * wk[k0] };
    // If maximum water capacity is not exceeded - update VolWaterContent of this cell and exit the function.
    if dripw[0] <= h2odef {
        unsafe {
            VolWaterContent[l0][k0] += dripw[0] / (dl[l0] * wk[k0]);
        }
        return Ok(());
    }
    // If maximum water capacity is exceeded - calculate the excess of water flowing out of this cell (in cm3 per slab) as dripw[1].
    // The next ring of cells (kr=1) will receive it as incoming water flow.
    dripw[1] = dripw[0] - h2odef;
    // Compute the movement of nitrate N to the next ring
    let mut cnw; // concentration of nitrate N in the outflowing water
    if unsafe { VolNo3NContent[l0][k0] } > 1.0e-30 {
        cnw = unsafe {
            VolNo3NContent[l0][k0] / (VolWaterContent[l0][k0] + dripw[0] / (dl[l0] * wk[k0]))
        };
        // cnw is multiplied by dripw[1] to get dripn[1], the amount of nitrate N going out to the next ring of cells.
        // It is assumed, however, that not more than a proportion (NO3FlowFraction) of the nitrate N in this cell can be removed in one iteration.
        if unsafe { (cnw * MaxWaterCapacity[l0]) < (NO3FlowFraction[l0] * VolNo3NContent[l0][k0]) }
        {
            dripn[1] = unsafe { NO3FlowFraction[l0] * VolNo3NContent[l0][k0] * dl[l0] * wk[k0] };
            unsafe { VolNo3NContent[l0][k0] = (1. - NO3FlowFraction[l0]) * VolNo3NContent[l0][k0] };
        } else {
            dripn[1] = dripw[1] * cnw;
            unsafe {
                VolNo3NContent[l0][k0] = MaxWaterCapacity[l0] * cnw;
            }
        }
    }
    // The movement of urea N to the next ring is computed similarly.
    let mut cuw; // concentration of urea N in the outflowing water
    if unsafe { VolUreaNContent[l0][k0] } > 1.0e-30 {
        cuw = unsafe {
            VolUreaNContent[l0][k0] / (VolWaterContent[l0][k0] + dripw[0] / (dl[l0] * wk[k0]))
        };
        if unsafe { (cuw * MaxWaterCapacity[l0]) < (NO3FlowFraction[l0] * VolUreaNContent[l0][k0]) }
        {
            dripu[1] = unsafe { NO3FlowFraction[l0] * VolUreaNContent[l0][k0] * dl[l0] * wk[k0] };
            unsafe {
                VolUreaNContent[l0][k0] = (1. - NO3FlowFraction[l0]) * VolUreaNContent[l0][k0];
            }
        } else {
            dripu[1] = dripw[1] * cuw;
            unsafe {
                VolUreaNContent[l0][k0] = MaxWaterCapacity[l0] * cuw;
            }
        }
    }
    let mut defcit = ndarray::Array2::<f64>::zeros((40, 20)); // array of the difference between water capacity and actual water content in each cell of the ring

    // Set VolWaterContent of the cell in which the drip is located to MaxWaterCapacity.
    unsafe {
        VolWaterContent[l0][k0] = MaxWaterCapacity[l0];
    }
    // Loop of concentric rings of cells, starting from ring 1.
    // Assign zero to the sums sv, st, sn, sn1, su and su1.
    for kr in 1..maxl as usize {
        let mut uplimit; //  upper limit of soil water content in a soil cell
        let mut sv = 0f64; // sum of actual water content in a ring of cells, cm3
        let mut st = 0f64; // sum of total water capacity in a ring of cells, cm3
        let mut sn = 0f64; // sum of nitrate N content in a ring of cells, mg.
        let mut sn1 = 0f64; // sum of movable nitrate N content in a ring of cells, mg
        let mut su = 0f64; // sum of urea N content in a ring of cells, mg
        let mut su1 = 0f64; // sum of movable urea N content in a ring of cells, mg
        let radius = (6 * kr) as f64; // radius (cm) of the wetting ring
        let mut dist; // distance (cm) of a cell center from drip location

        for l in 1..unsafe { nl as usize } {
            // Upper limit of water content is the porespace volume in layers below the water table, MaxWaterCapacity in other layers.
            if l >= unsafe { WaterTableLayer as usize } {
                uplimit = unsafe { PoreSpace[l] };
            } else {
                unsafe {
                    uplimit = MaxWaterCapacity[l];
                }
            }
            for k in 0..unsafe { nk as usize } {
                // Compute the sums sv, st, sn, sn1, su and su1 within the radius limits of this ring.
                // The array defcit is the sum of difference between uplimit and VolWaterContent of each cell.
                dist = utils::cell_distance(l, k, l0, k0, unsafe { RowSpace })?;
                if dist <= radius && dist > (radius - 6.) {
                    sv += unsafe { VolWaterContent[l][k] * dl[l] * wk[k] };
                    st += uplimit * unsafe { dl[l] * wk[k] };
                    sn += unsafe { VolNo3NContent[l][k] * dl[l] * wk[k] };
                    sn1 += unsafe { VolNo3NContent[l][k] * dl[l] * wk[k] * NO3FlowFraction[l] };
                    su += unsafe { VolUreaNContent[l][k] * dl[l] * wk[k] };
                    su1 += unsafe { VolUreaNContent[l][k] * dl[l] * wk[k] * NO3FlowFraction[l] };
                    defcit[[l, k]] = uplimit - unsafe { VolWaterContent[l][k] };
                } else {
                    defcit[[l, k]] = 0.;
                }
            }
        }
        // Compute the amount of water needed to saturate all the cells in this ring (h2odef).
        h2odef = st - sv;
        // Test if the amount of incoming flow, dripw(kr), is greater than h2odef.
        if dripw[kr] <= h2odef {
            // In this case, this will be the last wetted ring.
            // Update VolWaterContent in this ring, by wetting each cell in proportion to its defcit.
            // Update VolNo3NContent and VolUreaNContent of the cells in this ring by the same proportion.
            // This is executed for all the cells in the ring.
            for l in 1..unsafe { nl } as usize {
                for k in 0..unsafe { nk } as usize {
                    dist = utils::cell_distance(l, k, l0, k0, unsafe { RowSpace })?;
                    if dist <= radius && dist > (radius - 6.) {
                        unsafe {
                            VolWaterContent[l][k] += dripw[kr] * defcit[[l, k]] / h2odef;
                            VolNo3NContent[l][k] += dripn[kr] * defcit[[l, k]] / h2odef;
                            VolUreaNContent[l][k] += dripu[kr] * defcit[[l, k]] / h2odef;
                        }
                    }
                }
            }
            return Ok(());
        }
        // If dripw(kr) is greater than h2odef, calculate cnw and cuw as the concentration of nitrate and urea N in the total water of this ring after receiving the incoming water and nitrogen.
        cnw = (sn + dripn[kr]) / (sv + dripw[kr]);
        cuw = (su + dripu[kr]) / (sv + dripw[kr]);
        let drwout = dripw[kr] - h2odef; //  the amount of water going out of a ring, cm3.

        // Compute the nitrate and urea N going out of this ring, and their amount lost from this ring.
        // It is assumed that not more than a certain part of the total nitrate or urea N (previously computed as sn1 an su1) can be lost from a ring in one iteration.
        // drnout and xnloss are adjusted accordingly.
        // druout and xuloss are computed similarly for urea N.
        let mut drnout = drwout * cnw; //  the amount of nitrate N going out of a ring, mg
        let mut xnloss = 0f64; // the amount of nitrate N lost from a ring, mg
        if drnout > dripn[kr] {
            xnloss = drnout - dripn[kr];
            if xnloss > sn1 {
                xnloss = sn1;
                drnout = dripn[kr] + xnloss;
            }
        }
        let mut druout = drwout * cuw; //  the amount of urea N going out of a ring, mg
        let mut xuloss = 0f64; // the amount of urea N lost from a ring, mg
        if druout > dripu[kr] {
            xuloss = druout - dripu[kr];
            if xuloss > su1 {
                xuloss = su1;
                druout = dripu[kr] + xuloss;
            }
        }
        // For all the cells in the ring, as in the 1st cell, saturate VolWaterContent to uplimit, and update VolNo3NContent and VolUreaNContent.
        for l in 1..unsafe { nl } as usize {
            uplimit = unsafe {
                if l >= WaterTableLayer as usize {
                    PoreSpace[l]
                } else {
                    MaxWaterCapacity[l]
                }
            };
            for k in 0..unsafe { nk } as usize {
                dist = utils::cell_distance(l, k, l0, k0, unsafe { RowSpace })?;
                if dist <= radius && dist > (radius - 6.) {
                    unsafe {
                        VolWaterContent[l][k] = uplimit;
                        VolNo3NContent[l][k] = if xnloss <= 0. {
                            uplimit * cnw
                        } else {
                            VolNo3NContent[l][k] * (1. - xnloss / sn)
                        };
                        VolUreaNContent[l][k] = if xuloss <= 0. {
                            uplimit * cuw
                        } else {
                            VolUreaNContent[l][k] * (1. - xuloss / su)
                        };
                    }
                }
            }
        }
        // The outflow of water, nitrate and urea from this ring will be the inflow into the next ring.
        if kr < (unsafe { nl } as usize - l0 - 1) && kr < maxl as usize - 1 {
            dripw[kr + 1] = drwout;
            dripn[kr + 1] = drnout;
            dripu[kr + 1] = druout;
        } else
        // If this is the last ring, the outflowing water will be added to the drainage, CumWaterDrained, the outflowing nitrogen to SoilNitrogenLoss.
        {
            unsafe {
                CumWaterDrained += 10. * drwout / RowSpace;
                SoilNitrogenLoss += drnout + druout;
            }
            return Ok(());
        }
        // Repeat all these procedures for the next ring.
    }
    Ok(())
}
