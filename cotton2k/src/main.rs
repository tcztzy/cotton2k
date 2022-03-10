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
    let mut profile: cotton2k::Profile = toml::from_str(&contents)?;
    profile.path = profile_path.to_path_buf();
    match profile.weather_path.is_relative() {
        true => {
            profile.weather_path = profile_path.parent().unwrap().join(profile.weather_path);
        }
        false => {}
    }
    profile.soil_impedance = Some(profile_path.parent().unwrap().join("soil_imp.csv"));
    profile.initialize()?;
    profile.output_file_headers()?;
    unsafe {
        let count = (cotton2k::DayFinish - cotton2k::DayStart + 1) as u64;
        let mut pb = pbr::ProgressBar::new(count);
        // Do daily simulations
        cotton2k::Daynum = cotton2k::DayStart - 1;
        cotton2k::bEnd = false;
        // Start the daily loop. If variable bEnd has been assigned a value of true end simulation.
        while !cotton2k::bEnd {
            let need_to_adjust = profile.adjust();
            pb.inc();
            // Execute simulation for this day.
            profile.simulate_this_day();
            profile.write_record()?;
            // If there are pending plant adjustments, call WriteStateVariables() to write
            // state variables of this day in a scratch file.
            if need_to_adjust {
                cotton2k::WriteStateVariables(true);
            }
        }
    }
    Ok(())
}
