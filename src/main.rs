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
        for _ in cotton2k::DayStart..(cotton2k::DayFinish + 1) {
            let mut state = if profile.states.len() > 0 {
                let mut new_state = profile.states.last().unwrap().clone();
                new_state.date = new_state.date.succ();
                new_state
            } else {
                cotton2k::State::new(
                    &profile,
                    chrono::NaiveDate::from_yo(cotton2k::iyear, cotton2k::DayStart as u32),
                )
            };
            // let need_to_adjust = profile.adjust()?;
            pb.inc();
            // Execute simulation for this day.
            state.simulate_this_day(&mut profile)?;
            profile.write_record()?;
            // If there are pending plant adjustments, call WriteStateVariables() to write
            // state variables of this day in a scratch file.
            // if need_to_adjust {
            //     cotton2k::WriteStateVariables(true);
            // }

            if cotton2k::bEnd {
                break;
            }
        }
    }
    Ok(())
}
