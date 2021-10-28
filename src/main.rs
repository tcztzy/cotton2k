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
        write!(formatter, "a string represents chrono::NaiveDateTime")
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

fn from_isoformat<'de, D>(d: D) -> Result<NaiveDate, D::Error>
where
    D: de::Deserializer<'de>,
{
    d.deserialize_str(NaiveDateVisitor)
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
    #[serde(deserialize_with = "from_isoformat")]
    emerge_date: NaiveDate,
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
        DayEmerge = profile.emerge_date.ordinal() as i32;
        DayFinish = profile.stop_date.ordinal() as i32;
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
