#![feature(panic_always_abort)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
use chrono::{Datelike, NaiveDate};
use serde::{de, Deserialize};
use std::ffi::CString;
use std::io::Read;
use std::path::Path;

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

struct NaiveDateVisitor;

impl<'de> de::Visitor<'de> for NaiveDateVisitor {
    type Value = NaiveDate;

    fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(formatter, "a string represents chrono::NaiveDate")
    }

    fn visit_str<E>(self, s: &str) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        match NaiveDate::parse_from_str(s, "%F") {
            Ok(t) => Ok(t),
            Err(_) => Err(de::Error::invalid_value(de::Unexpected::Str(s), &self)),
        }
    }
}
struct OptionalNaiveDateVisitor;

impl<'de> de::Visitor<'de> for OptionalNaiveDateVisitor {
    type Value = Option<NaiveDate>;

    fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(formatter, "a optional string represents chrono::NaiveDate")
    }

    fn visit_none<E>(self) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        Ok(None)
    }

    fn visit_some<D>(self, deserializer: D) -> Result<Self::Value, D::Error>
    where
        D: de::Deserializer<'de>,
    {
        Ok(Some(deserializer.deserialize_str(NaiveDateVisitor)?))
    }
}

fn from_isoformat<'de, D>(d: D) -> Result<NaiveDate, D::Error>
where
    D: de::Deserializer<'de>,
{
    d.deserialize_str(NaiveDateVisitor)
}

fn from_isoformat_option<'de, D>(d: D) -> Result<Option<NaiveDate>, D::Error>
where
    D: de::Deserializer<'de>,
{
    d.deserialize_option(OptionalNaiveDateVisitor)
}

#[derive(Deserialize, Debug)]
struct Profile {
    name: Option<String>,
    latitude: f64,
    longitude: f64,
    elevation: f64,
    #[serde(deserialize_with = "from_isoformat")]
    start_date: NaiveDate,
    #[serde(deserialize_with = "from_isoformat")]
    stop_date: NaiveDate,
    #[serde(default)]
    #[serde(deserialize_with = "from_isoformat_option")]
    emerge_date: Option<NaiveDate>,
    #[serde(default)]
    #[serde(deserialize_with = "from_isoformat_option")]
    plant_date: Option<NaiveDate>,
    co2_enrichment: Option<CO2Enrichment>,
    mulch: Option<Mulch>,
}

#[derive(Deserialize, Debug)]
struct CO2Enrichment {
    factor: f64,
    #[serde(deserialize_with = "from_isoformat")]
    start_date: NaiveDate,
    #[serde(deserialize_with = "from_isoformat")]
    stop_date: NaiveDate,
}

#[derive(Deserialize, Debug)]
enum MulchType {
    NoMulch,
    All,                  // plastic layer on all soil surface
    OneColumnLeftAtSide, // plastic layer on all soil surface except one column at each side of the plant row.
    TwoColumnsLeftAtSide, // plastic layer on all soil surface except two columns at each side of the plant row.
}

#[derive(Deserialize, Debug)]
struct Mulch {
    indicator: MulchType,
    sw_trans: f64,
    lw_trans: f64,
    #[serde(deserialize_with = "from_isoformat")]
    start_date: NaiveDate,
    #[serde(default)]
    #[serde(deserialize_with = "from_isoformat_option")]
    stop_date: Option<NaiveDate>,
}

fn read_profile(profile_path: &Path) {
    let mut file = std::fs::File::open(profile_path).unwrap();
    let mut contents = String::new();
    file.read_to_string(&mut contents).unwrap();
    let mut profile: Profile = toml::from_str(&contents).unwrap();
    match profile.name {
        None => {
            profile.name = Some(String::from(
                profile_path.file_stem().unwrap().to_str().unwrap(),
            ));
        }
        Some(_) => {}
    }
    unsafe {
        Latitude = profile.latitude;
        Longitude = profile.longitude;
        Elevation = profile.elevation;
        iyear = profile.start_date.year();
        DayStart = profile.start_date.ordinal() as i32;
        DayFinish = profile.stop_date.ordinal() as i32;
        match profile.emerge_date {
            Some(date) => DayEmerge = date.ordinal() as i32,
            None => {}
        }
        match profile.plant_date {
            Some(date) => {
                DayPlant = date.ordinal() as i32;
            }
            None => {}
        }
        if profile.emerge_date.is_none() {
            // If the date of emergence has not been given, emergence will be simulated by the model. In this case, isw = 0, and a check is performed to make sure that the date of planting has been given.
            if profile.plant_date.is_none() {
                panic!(
                    "one of planting date or emergence date must be given in the profile file!!"
                );
            }
            isw = 0;
        } else if profile.emerge_date > profile.plant_date {
            // If the date of emergence has been given in the input: isw = 1 if simulation starts before emergence,
            isw = 1;
        } else {
            // or isw = 2 if simulation starts at emergence.
            isw = 2;
            Kday = 1;
        }
        // For advanced users only: if there is CO2 enrichment, read also CO2 factor, DOY dates
        match profile.co2_enrichment {
            Some(enrichment) => {
                CO2EnrichmentFactor = enrichment.factor;
                DayStartCO2 = enrichment.start_date.ordinal() as i32;
                DayEndCO2 = enrichment.stop_date.ordinal() as i32;
            }
            None => CO2EnrichmentFactor = 0.,
        }
        // If soil mulch is used, read relevant parameters.
        match profile.mulch {
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
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        panic!("profile file path should be provided!");
    }
    let mut app: C2KApp = unsafe { C2KApp::new() };
    read_profile(Path::new(&args[1]));
    let filename = CString::new(args[2].to_string()).expect("Error");
    unsafe {
        app.RunTheModel(filename.as_ptr());
    }
}
