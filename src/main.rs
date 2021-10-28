#![feature(panic_always_abort)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
use chrono::{Datelike, NaiveDate};
use serde::{de, Deserialize};
use std::ffi::CString;
use std::io::Read;
use std::path::{Path, PathBuf};

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
    weather_path: PathBuf,
    site: Site,
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

#[derive(Deserialize, Debug)]
struct Site {
    average_wind_speed: Option<f64>,
    estimate_dew_point: (f64, f64),
}

#[derive(Deserialize, Debug)]
struct WeatherRecord {
    #[serde(deserialize_with = "from_isoformat")]
    date: NaiveDate,
    irradiation: f64,
    tmax: f64,
    tmin: f64,
    rain: f64,
    wind: Option<f64>,
    tdew: Option<f64>,
}

fn read_profile(profile_path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = std::fs::File::open(profile_path)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    let mut profile: Profile = toml::from_str(&contents)?;
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
    let weather_path = if profile.weather_path.is_relative() {
        profile_path.parent().unwrap().join(profile.weather_path)
    } else {
        profile.weather_path
    };
    println!("{:?}", weather_path.to_str().unwrap());
    let mut rdr = csv::Reader::from_path(weather_path)?;
    let mut jdd: i32 = 0;
    for result in rdr.deserialize() {
        let record: WeatherRecord = result?;
        jdd = record.date.ordinal() as i32;
        let j = jdd - unsafe { DayStart };
        if j < 0 {continue;}
        unsafe {
            Clim[j as usize].nDay = jdd;
            // convert \frac{MJ}{m^2} to langleys
            Clim[j as usize].Rad = record.irradiation * 23.884;
            Clim[j as usize].Tmax = record.tmax;
            Clim[j as usize].Tmin = record.tmin;
            Clim[j as usize].Wind = if profile.site.average_wind_speed.is_some() && record.wind.is_none() {
                profile.site.average_wind_speed.unwrap()
            } else {
                record.wind.unwrap_or(0.)
            };
            Clim[j as usize].Tdew = record.tdew.unwrap_or(estimate_dew_point(
                record.tmax,
                profile.site.estimate_dew_point.0,
                profile.site.estimate_dew_point.1,
            ));
            Clim[j as usize].Rain = record.rain;
        }
    }
    unsafe {
        LastDayWeatherData = jdd;
    }
    Ok(())
}

/// estimates the approximate daily average dewpoint temperature when it is not available.
pub fn estimate_dew_point(maxt: f64, site5: f64, site6: f64) -> f64 {
    if maxt <= 20. {
        site5
    } else if maxt >= 40. {
        site6
    } else {
        ((40. - maxt) * site5 + (maxt - 20.) * site6) / 20.
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        panic!("profile file path should be provided!");
    }
    let mut app: C2KApp = unsafe { C2KApp::new() };
    read_profile(Path::new(&args[1])).expect("Error in read_profile");
    let filename = CString::new(args[2].to_string()).expect("Error");
    unsafe {
        app.RunTheModel(filename.as_ptr());
    }
}
