use crate::bindings::*;
use crate::profile::{
    AgronomyOperation, FertilizationMethod, IrrigationMethod, MulchType, Profile, WeatherRecord,
};
use chrono::Datelike;
use std::io::Read;
use std::path::Path;

pub fn read_profile(profile_path: &Path) -> Result<(), Box<dyn std::error::Error>> {
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
        InitializeGlobal();
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
        SitePar[1] = profile.site.wind_blow_after_sunrise;
        SitePar[2] = profile.site.wind_max_after_noon;
        SitePar[3] = profile.site.wind_stop_after_sunset;
        SitePar[4] = profile.site.night_time_wind_factor;
        SitePar[7] = profile.site.cloud_type_correction_factor;
        SitePar[8] = profile.site.max_temperature_after_noon;
        SitePar[9] = profile.site.deep_soil_temperature.0;
        SitePar[10] = profile.site.deep_soil_temperature.1;
        SitePar[11] = profile.site.deep_soil_temperature.2;
        SitePar[12] = profile.site.dew_point_range.0;
        SitePar[13] = profile.site.dew_point_range.1;
        SitePar[14] = profile.site.dew_point_range.2;
        SitePar[15] = profile.site.albedo_range.1;
        SitePar[16] = profile.site.albedo_range.0;
        for pair in profile.cultivar_parameters.iter().enumerate() {
            VarPar[pair.0 + 1] = *pair.1;
        }
        RowSpace = if profile.skip_row_width > 0. {
            (profile.row_space + profile.skip_row_width) / 2.
        } else {
            profile.row_space
        };
        // PlantRowLocation is the distance from edge of slab, cm, of the plant row.
        PlantRowLocation = RowSpace / 2.;
        // Compute PlantPopulation - number of plants per hectar, and
        // PerPlantArea - the average surface area per plant, in dm2, and
        // the empirical plant density factor (DensityFactor). This factor will be used to express the effect of plant density on some plant growth rate functions.  Note that DensityFactor =1 for 5 plants per sq m (or 50000 per ha).
        PlantPopulation = profile.plants_per_meter / RowSpace * 1000000.;
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
    let weather_path = if profile.weather_path.is_relative() {
        profile_path.parent().unwrap().join(profile.weather_path)
    } else {
        profile.weather_path
    };
    let mut rdr = csv::Reader::from_path(weather_path)?;
    let mut jdd: i32 = 0;
    for result in rdr.deserialize() {
        let record: WeatherRecord = result?;
        jdd = record.date.ordinal() as i32;
        let j = jdd - unsafe { DayStart };
        if j < 0 {
            continue;
        }
        unsafe {
            Clim[j as usize].nDay = jdd;
            // convert \frac{MJ}{m^2} to langleys
            Clim[j as usize].Rad = record.irradiation * 23.884;
            Clim[j as usize].Tmax = record.tmax;
            Clim[j as usize].Tmin = record.tmin;
            Clim[j as usize].Wind =
                if profile.site.average_wind_speed.is_some() && record.wind.is_none() {
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
    let mut idef: usize = 0;
    let mut icult: usize = 0;
    let mut ipix: usize = 0;
    unsafe {
        NumNitApps = 0;
        NumIrrigations = 0;
        NumWaterTableData = 0;
        for i in 0..5 {
            DefoliationDate[i] = 0;
            DefoliationMethod[i] = 0;
            DefoliantAppRate[i] = 0.;
        }

        for ao in profile.agronomy_operations {
            match ao {
                AgronomyOperation::irrigation {
                    date,
                    amount,
                    predict,
                    method,
                    drip_x,
                    drip_y,
                    max_amount,
                    stop_predict_date,
                } => {
                    if predict {
                        MaxIrrigation = max_amount.unwrap();
                        DayStartPredIrrig = date.ordinal() as i32;
                        DayStopPredIrrig = stop_predict_date.unwrap().ordinal() as i32;
                        if let IrrigationMethod::Drip = method {
                            LocationColumnDrip = slab_horizontal_location(drip_x) as i32;
                            LocationLayerDrip = slab_vertical_location(drip_y) as i32;
                        }
                        IrrigMethod = method as i32;
                    } else {
                        Irrig[NumIrrigations as usize].day = date.ordinal() as i32;
                        Irrig[NumIrrigations as usize].amount = amount;
                        if let IrrigationMethod::Drip = method {
                            Irrig[NumIrrigations as usize].LocationColumnDrip =
                                slab_horizontal_location(drip_x) as i32;
                            Irrig[NumIrrigations as usize].LocationLayerDrip =
                                slab_vertical_location(drip_y) as i32;
                        }
                        Irrig[NumIrrigations as usize].method = method as i32;
                        NumIrrigations += 1;
                    }
                }
                AgronomyOperation::fertilization {
                    date,
                    urea,
                    nitrate,
                    ammonium,
                    method,
                    drip_x,
                    drip_y,
                } => {
                    NFertilizer[NumNitApps as usize].day = date.ordinal() as i32;
                    NFertilizer[NumNitApps as usize].amtamm = ammonium;
                    NFertilizer[NumNitApps as usize].amtnit = nitrate;
                    NFertilizer[NumNitApps as usize].amtura = urea;
                    match method {
                        FertilizationMethod::Sidedress | FertilizationMethod::Drip => {
                            NFertilizer[NumNitApps as usize].ksdr =
                                slab_horizontal_location(drip_x) as i32;
                            NFertilizer[NumNitApps as usize].lsdr =
                                slab_vertical_location(drip_y) as i32;
                        }
                        _ => {}
                    }
                    NFertilizer[NumNitApps as usize].mthfrt = method as i32;
                    NumNitApps += 1;
                }
                AgronomyOperation::defoliation {
                    date,
                    open_ratio,
                    predict,
                    ppa,
                } => {
                    DefoliationDate[idef] = date.ordinal() as i32;
                    DefoliantAppRate[idef] = if predict { -99.9 } else { ppa };
                    DefoliationMethod[idef] = open_ratio;
                    DayFirstDef = DefoliationDate[0];
                    idef += 1;
                }
                AgronomyOperation::cultivation { date, depth } => {
                    CultivationDate[icult] = date.ordinal() as i32;
                    CultivationDepth[icult] = depth;
                    icult += 1;
                }
                AgronomyOperation::pix { date, method, ppa } => {
                    pixday[ipix] = date.ordinal() as i32;
                    pixmth[ipix] = method as i32;
                    pixppa[ipix] = ppa;
                    ipix += 1;
                }
                AgronomyOperation::watertable { date, level, ecs } => {
                    DayWaterTableInput[NumWaterTableData as usize] = date.ordinal() as i32;
                    LevelsOfWaterTable[NumWaterTableData as usize] = level;
                    ElCondSatSoil[NumWaterTableData as usize] = ecs;
                    NumWaterTableData += 1;
                }
            }
        }
    }
    Ok(())
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

fn slab_vertical_location(distance: f64) -> usize {
    let mut sumdl = 0.;
    let mut result: usize = 0;
    unsafe {
        for l in 0..nl as usize {
            sumdl += dl[l];
            if sumdl >= distance {
                result = l;
                break;
            }
        }
    }
    result
}

fn slab_horizontal_location(distance: f64) -> usize {
    let mut sumwk = 0.;
    let mut result: usize = 0;
    unsafe {
        for k in 0..nk as usize {
            sumwk += wk[k];
            if sumwk >= distance {
                result = k;
                break;
            }
        }
    }
    result
}
unsafe fn InitializeGlobal()
//     This function initializes many "global" variables at the start of a
//  simulation. It is called from ReadInput(). Note that initialization
//  is needed at the start of each simulation (NOT at start of the run).
{
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
    CO2EnrichmentFactor = 0.;
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
    NumWaterTableData = 0;
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
    TotalLeafArea = 0.;
    TotalLeafWeight = 0.20;
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
    }
    //
    for i in 0..20 {
        AbscissionLag[i] = 0.;
        DayWaterTableInput[i] = 0;
        ElCondSatSoil[i] = 0.;
        LevelsOfWaterTable[i] = 0.;
        ShedByCarbonStress[i] = 0.;
        ShedByNitrogenStress[i] = 0.;
        ShedByWaterStress[i] = 0.;
    }
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
