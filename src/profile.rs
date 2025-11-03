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
                let mut state = if self.states.len() > 0 {
                    let mut new_state = self.states.last().unwrap().clone();
                    new_state.date = new_state.date.succ_opt().unwrap();
                    new_state
                } else {
                    State::new(
                        self,
                        NaiveDate::from_yo_opt(iyear, DayStart as u32).unwrap(),
                    )
                };
                // Execute simulation for this day.
                match state.simulate_this_day(self) {
                    Err(e) => {
                        if e.level == 0 {
                            println!("{}", e.message);
                            break;
                        }
                    }
                    _ => {}
                }
                self.write_record()?;
                self.states.push(state);
                if bEnd {
                    break;
                }
            }
        }
        Ok(())
    }

    pub fn initialize(self: &mut Self) -> Result<(), Box<dyn std::error::Error>> {
        self.lmax = 0.;
        self.num_watertable_data = 0;
        unsafe {
            InitializeGlobal();
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
        let mut jdd: u32 = 0;
        for result in rdr.deserialize() {
            let record: WeatherRecord = result?;
            jdd = record.date.ordinal();
            let j = jdd as i32 - unsafe { DayStart };
            if j < 0 {
                continue;
            }
            unsafe {
                Clim[j as usize].nDay = jdd as i32;
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
        self.last_day_weather_data = NaiveDate::from_yo_opt(unsafe { iyear }, jdd).unwrap();
        let mut idef: usize = 0;
        let mut icult: usize = 0;
        unsafe {
            NumIrrigations = 0;
            for i in 0..5 {
                DefoliationDate[i] = 0;
                DefoliationMethod[i] = 0;
                DefoliantAppRate[i] = 0.;
            }

            for ao in &self.agronomy_operations {
                match ao {
                    AgronomyOperation::irrigation(irrigation) => {
                        if irrigation.predict {
                            MaxIrrigation = irrigation.max_amount.unwrap();
                            DayStartPredIrrig = irrigation.date.ordinal() as i32;
                            DayStopPredIrrig =
                                irrigation.stop_predict_date.unwrap().ordinal() as i32;
                            if let IrrigationMethod::Drip = irrigation.method {
                                LocationColumnDrip =
                                    utils::slab_horizontal_location(irrigation.drip_x, RowSpace)?
                                        as i32;
                                LocationLayerDrip =
                                    utils::slab_vertical_location(irrigation.drip_y)? as i32;
                            }
                            IrrigMethod = irrigation.method as i32;
                        } else {
                            Irrig[NumIrrigations as usize].day = irrigation.date.ordinal() as i32;
                            Irrig[NumIrrigations as usize].amount = irrigation.amount;
                            if let IrrigationMethod::Drip = irrigation.method {
                                Irrig[NumIrrigations as usize].LocationColumnDrip =
                                    utils::slab_horizontal_location(irrigation.drip_x, RowSpace)?
                                        as i32;
                                Irrig[NumIrrigations as usize].LocationLayerDrip =
                                    utils::slab_vertical_location(irrigation.drip_y)? as i32;
                            }
                            Irrig[NumIrrigations as usize].method = irrigation.method as i32;
                            NumIrrigations += 1;
                        }
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
                    _ => {}
                }
            }
            self.num_watertable_data = self
                .agronomy_operations
                .iter()
                .filter(|x| match x {
                    AgronomyOperation::watertable { .. } => true,
                    _ => false,
                })
                .count();
            if matches!(self.light_intercept_method, LightInterceptMethod::Latered)
                && self.light_intercept_parameters.is_none()
            {
                panic!(
                    "light_intercept_parameters must be provided when using the Latered light intercept method"
                );
            }
        }
        self.read_soil_impedance(self.soil_impedance.as_ref().unwrap())?;
        unsafe {
            InitSoil(&self.soil_layers, &self.soil_hydraulic);
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

    pub fn write_record(self: &Self) -> Result<(), Box<dyn std::error::Error>> {
        let mut f = std::fs::OpenOptions::new()
            .write(true)
            .append(true)
            .open(self.path.parent().unwrap().join("output.csv"))?;
        let mut record = vec![
            chrono::NaiveDate::from_yo_opt(unsafe { iyear }, unsafe { Daynum } as u32)
                .unwrap()
                .format("%F")
                .to_string(),
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
unsafe fn InitializeRootData() {
    // The parameters of the root model are defined for each root class:
    // grind(i), cuind(i), thtrn(i), trn(i), thdth(i), dth(i).
    NumRootAgeGroups = 3;
    cgind[0] = 1.;
    cgind[1] = 1.;
    cgind[2] = 0.10;
    let rlint = 10.; // Vertical interval, in cm, along the taproot, for initiating lateral roots.
    let mut ll = 1; // Counter for layers with lateral roots.
    let mut sumdl = 0.; // Distance from soil surface to the middle of a soil layer.
    for l in 0..nl as usize {
        // Using the value of rlint (interval between lateral roots), the layers from which lateral roots may be initiated are now computed.
        // LateralRootFlag[l] is assigned a value of 1 for these layers.
        LateralRootFlag[l] = 0;
        if l > 0 {
            sumdl += 0.5 * dl[l - 1];
        }
        sumdl += 0.5 * dl[l];
        if sumdl >= ll as f64 * rlint {
            LateralRootFlag[l] = 1;
            ll += 1;
        }
    }
    // All the state variables of the root system are initialized to zero.
    for l in 0..nl as usize {
        if l < 3 {
            RootColNumLeft[l] = PlantRowColumn - 1;
            RootColNumRight[l] = PlantRowColumn + 2;
        } else if l < 7 {
            RootColNumLeft[l] = PlantRowColumn;
            RootColNumRight[l] = PlantRowColumn + 1;
        } else {
            RootColNumLeft[l] = 0;
            RootColNumRight[l] = 0;
        }
        //
        for k in 0..nk as usize {
            PotGroRoots[l][k] = 0.;
            RootGroFactor[l][k] = 1.;
            ActualRootGrowth[l][k] = 0.;
            RootAge[l][k] = 0.;
            for i in 0..3 {
                RootWeight[l][k][i] = 0.;
            }
        }
    }
    //
    RootWeight[0][(PlantRowColumn - 1) as usize][0] = 0.0020;
    RootWeight[0][PlantRowColumn as usize][0] = 0.0070;
    RootWeight[0][PlantRowColumn as usize + 1][0] = 0.0070;
    RootWeight[0][PlantRowColumn as usize + 2][0] = 0.0020;
    RootWeight[1][(PlantRowColumn - 1) as usize][0] = 0.0040;
    RootWeight[1][PlantRowColumn as usize][0] = 0.0140;
    RootWeight[1][PlantRowColumn as usize + 1][0] = 0.0140;
    RootWeight[1][PlantRowColumn as usize + 2][0] = 0.0040;
    RootWeight[2][(PlantRowColumn - 1) as usize][0] = 0.0060;
    RootWeight[2][PlantRowColumn as usize][0] = 0.0210;
    RootWeight[2][PlantRowColumn as usize + 1][0] = 0.0210;
    RootWeight[2][PlantRowColumn as usize + 2][0] = 0.0060;
    RootWeight[3][PlantRowColumn as usize][0] = 0.0200;
    RootWeight[3][PlantRowColumn as usize + 1][0] = 0.0200;
    RootWeight[4][PlantRowColumn as usize][0] = 0.0150;
    RootWeight[4][PlantRowColumn as usize + 1][0] = 0.0150;
    RootWeight[5][PlantRowColumn as usize][0] = 0.0100;
    RootWeight[5][PlantRowColumn as usize + 1][0] = 0.0100;
    RootWeight[6][PlantRowColumn as usize][0] = 0.0050;
    RootWeight[6][PlantRowColumn as usize + 1][0] = 0.0050;
    // Start loop for all soil layers containing roots.
    DepthLastRootLayer = 0.;
    TotalRootWeight = 0.;
    for l in 0..7 {
        DepthLastRootLayer += dl[l]; // compute total depth to the last layer with roots (DepthLastRootLayer).
        for k in 0..nk as usize {
            // For each soil soil cell with roots, compute total root weight per plant (TotalRootWeight), and convert RootWeight from g per plant to g per cell.
            for i in 0..3 {
                TotalRootWeight += RootWeight[l][k][i];
                RootWeight[l][k][i] = RootWeight[l][k][i] * 0.01 * RowSpace / PerPlantArea;
            }
            // initialize RootAge to a non-zero value for each cell containing roots.
            if RootWeight[l][k][0] > 0. {
                RootAge[l][k] = 0.01;
            }
        }
    }
    //     Initial value of taproot length, TapRootLength, is computed to the
    // middle of the last layer with roots. The last soil layer with
    // taproot, LastTaprootLayer, is defined.
    NumLayersWithRoots = 7;
    TapRootLength = DepthLastRootLayer - 0.5 * dl[(NumLayersWithRoots - 1) as usize];
    LastTaprootLayer = 6;
}
