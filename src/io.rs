use super::util::DoyToDate;
use super::{Simulation, State};
use chrono::prelude::*;
use serde::ser::{Serialize, SerializeStruct, Serializer};
use std::ffi::CStr;
use std::fs::{File, OpenOptions};
use std::io::prelude::*;
use std::os::raw::c_char;

#[no_mangle]
extern "C" fn b01(profile_name: *const c_char, description: *const c_char) {
    let local: DateTime<Local> = Local::now();
    let profile_name = unsafe { CStr::from_ptr(profile_name).to_str().unwrap() };
    let description = unsafe { CStr::from_ptr(description).to_str().unwrap() };
    let mut file = File::create(format!("output/{}.B01", profile_name)).unwrap();
    writeln!(file, "{:>50}", "COTTON2K Version 4.0 (2003)").unwrap();
    writeln!(
        file,
        "{:>62}",
        "A simulation model for irrigated cotton in arid regions"
    )
    .unwrap();
    writeln!(file, "{:>50}", "Written by Avishalom Marani").unwrap();
    writeln!(file).unwrap();
    writeln!(file, "Profile Name:    {:<20}", profile_name).unwrap();
    writeln!(
        file,
        "Simulation Date: {:<30}",
        local.format("%A, %B %d, %Y").to_string()
    )
    .unwrap();
    writeln!(file, "Description:     {:<55}", description).unwrap();
    writeln!(file).unwrap();
}

/// This function writes the input data to File20 (*.B01). It is executed once
/// at the beginning of the simulation. It is called by ReadInput().
#[no_mangle]
extern "C" fn WriteInitialInputData(
    sim: &Simulation,
    uscs: bool,
    plants_per_meter: f64,
    skip_row_width: f64,
    plant_population: f64,
    actual_weather_file_name: *const c_char,
    last_day_of_actual_weather: i32,
    predicted_weather_file_name: *const c_char,
    cultural_input_file_name: *const c_char,
    initial_soil_data_file_name: *const c_char,
    soil_hydrology_file_name: *const c_char,
    site_name: *const c_char,
    var_name: *const c_char,
) {
    let profile_name = unsafe { CStr::from_ptr(sim.profile_name).to_str().unwrap() };
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(format!("output/{}.B01", profile_name))
        .unwrap();
    write!(
        file,
        "    Latitude.. {:>8.2}{:>7}",
        sim.latitude,
        if sim.latitude < 0f64 {
            "South"
        } else {
            "North"
        }
    )
    .unwrap();
    writeln!(
        file,
        "         Longitude.. {:>8.2}{:>7}",
        sim.longitude,
        if sim.longitude < 0f64 { "West" } else { "East" }
    )
    .unwrap();
    writeln!(
        file,
        "    Elevation ({}).... {:>8.2}",
        if uscs { "ft" } else { "m" },
        sim.elevation * if uscs { 3.28 } else { 1.0 }
    )
    .unwrap();
    writeln!(
        file,
        "    Start Simulation {}       Stop Simulation.... {}",
        unsafe {
            CStr::from_ptr(DoyToDate(sim.day_start as i32, sim.year) as *const i8)
                .to_str()
                .unwrap()
        },
        unsafe {
            CStr::from_ptr(DoyToDate(sim.day_finish as i32, sim.year) as *const i8)
                .to_str()
                .unwrap()
        }
    )
    .unwrap();
    write!(file, "    Planting date....{}", unsafe {
        CStr::from_ptr(DoyToDate(sim.day_plant as i32, sim.year) as *const i8)
            .to_str()
            .unwrap()
    })
    .unwrap();
    if sim.day_emerge == 0 {
        writeln!(file, "       Emergence date is simulated   ").unwrap();
    } else {
        writeln!(file, "       Emergence date...   {}", unsafe {
            CStr::from_ptr(DoyToDate(sim.day_emerge as i32, sim.year) as *const i8)
                .to_str()
                .unwrap()
        })
        .unwrap()
    }

    if !uscs {
        writeln!(
            file,
            "    Row Spacing (cm) {:>8.2}          Plants Per Row-m... {:>8.2}",
            sim.row_space, plants_per_meter
        )
        .unwrap();
        writeln!(
            file,
            "    Skip Width (cm). {:>8.2}          Plants Per Ha...... {:>8.1}",
            skip_row_width, plant_population
        )
        .unwrap();
    } else {
        writeln!(
            file,
            "    Row Spacing (in) {:>8.2}          Plants Per Row-ft.. {:>8.2}",
            sim.row_space / 2.54,
            plants_per_meter * 0.305
        )
        .unwrap();
        writeln!(
            file,
            "    Skip Width (in). {:>8.2}          Plants Per Acre.... {:>8.1}",
            skip_row_width / 2.54,
            plant_population * 0.405
        )
        .unwrap();
    };
    if sim.co2_enrichment_factor > 1.0 {
        writeln!(
            file,
            "          CO2 enrichment factor              ...  {:>8.4}",
            sim.co2_enrichment_factor
        )
        .unwrap();
        writeln!(
            file,
            "    from ...... {}    to ...... {}",
            unsafe {
                CStr::from_ptr(DoyToDate(sim.day_start_co2 as i32, sim.year) as *const i8)
                    .to_str()
                    .unwrap()
            },
            unsafe {
                CStr::from_ptr(DoyToDate(sim.day_end_co2 as i32, sim.year) as *const i8)
                    .to_str()
                    .unwrap()
            }
        )
        .unwrap();
    }
    /*
    if (m_mulchdata.length() > 0 && MulchIndicator > 0)
    {
        File20 << "   Polyethylene mulch cover. Transmissivity values are: " << endl;
        File20 << " For short waves:  ";
        File20.width(8);
        File20.precision(3);
        File20 << MulchTranSW;
        File20 << " For long waves:  ";
        File20.width(8);
        File20.precision(3);
        File20 << MulchTranLW << endl;
        File20 << " From Day of Year  ";
        File20.width(4);
        File20 << DayStartMulch;
        File20 << " to Day of Year  ";
        File20.width(4);
        File20 << DayEndMulch << endl;
        if (MulchIndicator == 1)
            File20 << " All soil surface covered by mulch." << endl;
        else
        {
            File20.width(6);
            File20.precision(2);
            if (MulchIndicator == 2)
                File20 << RowSpace / maxk;
            else if (MulchIndicator == 3)
                File20 << 2 * RowSpace / maxk;
            File20 << " cm on each side of planr rows not covered by mulch." << endl;
        }
    }*/
    let actual_weather_file_name =
        unsafe { CStr::from_ptr(actual_weather_file_name).to_str().unwrap() };
    if !actual_weather_file_name.is_empty() {
        writeln!(
            file,
            "    Actual Weather Input File:     {}",
            actual_weather_file_name
        )
        .unwrap();
        writeln!(
            file,
            "    Last date read from Actual Weather File: {}",
            unsafe {
                CStr::from_ptr(DoyToDate(last_day_of_actual_weather, sim.year) as *const i8)
                    .to_str()
                    .unwrap()
            }
        )
        .unwrap();
    }
    let predicted_weather_file_name = unsafe {
        CStr::from_ptr(predicted_weather_file_name)
            .to_str()
            .unwrap()
    };
    if !predicted_weather_file_name.is_empty() {
        writeln!(
            file,
            "    Predicted Weather Input File:  {}",
            predicted_weather_file_name
        )
        .unwrap();
    }
    writeln!(file, "    Cultural Input File:           {}", unsafe {
        CStr::from_ptr(cultural_input_file_name).to_str().unwrap()
    })
    .unwrap();
    writeln!(file, "    Initial Soil Data Input File:  {}", unsafe {
        CStr::from_ptr(initial_soil_data_file_name)
            .to_str()
            .unwrap()
    })
    .unwrap();
    writeln!(file, "    Soil Hydrology Input File:     {}", unsafe {
        CStr::from_ptr(soil_hydrology_file_name).to_str().unwrap()
    })
    .unwrap();
    writeln!(file, "    Site...     {}", unsafe {
        CStr::from_ptr(site_name).to_str().unwrap()
    })
    .unwrap();
    writeln!(file, "    Variety...  {}", unsafe {
        CStr::from_ptr(var_name).to_str().unwrap()
    })
    .unwrap();
    writeln!(file).unwrap();
}

#[no_mangle]
extern "C" fn b01_append_defoliant(profile_name: *const c_char, date: *const c_char) {
    let profile_name = unsafe { CStr::from_ptr(profile_name).to_str().unwrap() };
    let date = unsafe { CStr::from_ptr(date).to_str().unwrap() };
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(format!("output/{}.B01", profile_name))
        .unwrap();
    writeln!(file, " ****   defoliant applied on {}    ****", date).unwrap();
}
