#[pyo3::pymodule]
mod cotton2k {
    use pyo3::prelude::*;
    use std::io::Read;

    #[pyfunction]
    fn run(path: &str) -> PyResult<()> {
        let profile_path = std::path::Path::new(path);
        let mut file = std::fs::File::open(profile_path)?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)?;
        let mut profile: cotton2k::Profile = match profile_path
            .extension()
            .unwrap()
            .to_str()
            .unwrap()
            .to_lowercase()
            .as_str()
        {
            "toml" => toml::from_str(&contents).expect("Parse file error!"),
            _ => serde_json::from_str(&contents).expect("Parse file error!"),
        };
        profile.path = profile_path.to_path_buf();
        match profile.weather_path.is_relative() {
            true => {
                profile.weather_path = profile_path.parent().unwrap().join(profile.weather_path);
            }
            false => {}
        }
        profile.soil_impedance = Some(profile_path.parent().unwrap().join("soil_imp.csv"));
        profile.run().expect("Simulation error!");
        Ok(())
    }
}
