use crate::bindings::*;
use std::io::Write;

pub fn output_file_headers() -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = csv::Writer::from_path("output.csv")?;
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

pub fn write_record() -> Result<(), Box<dyn std::error::Error>> {
    let mut f = std::fs::OpenOptions::new()
        .write(true)
        .append(true)
        .open("output.csv")?;
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
            ((CottonWeightOpenBolls + CottonWeightGreenBolls) * PlantPopulation / 1000.).to_string()
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
