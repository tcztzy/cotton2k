use cotton2k::{profile::Profile, run};
use std::io::Read;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        panic!("profile file path should be provided!");
    }
    let profile_path = std::path::Path::new(&args[1]);
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
    match profile.weather_path.is_relative() {
        true => {
            profile.weather_path = profile_path.parent().unwrap().join(profile.weather_path);
        }
        false => {}
    }
    profile.soil_impedance = Some(profile_path.parent().unwrap().join("soil_imp.csv"));
    run(profile)
}
