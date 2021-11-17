use crate::bindings::*;
use crate::profile::{
    AgronomyOperation, FertilizationMethod, IrrigationMethod, MulchType, PlantMap, Profile,
    SoilHydraulic, SoilLayer, WeatherRecord,
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
        light_intercept_method = profile.light_intercept_method;
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
        if profile.light_intercept_method == 2 {
            light_intercept_parameter = 0.;
            for i in 0..20 {
                light_intercept_parameters[i] = profile.light_intercept_parameters.unwrap()[i];
                light_intercept_parameter += light_intercept_parameters[i];
            }
        }
    }
    read_soil_impedance(&profile_path.parent().unwrap().join("soil_imp.csv"))?;
    unsafe {
        InitSoil(profile.soil_layers, profile.soil_hydraulic);
        if profile.plant_maps.is_some() {
            ReadPlantMapInput(profile.plant_maps.unwrap());
        }
        InitializeRootData();
        //     initialize some variables at the start of simulation.
        SoilNitrogenAtStart = TotalSoilNo3N + TotalSoilNh4N + TotalSoilUreaN;
        PlantWeightAtStart = TotalRootWeight + TotalStemWeight + TotalLeafWeight() + ReserveC;
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
            NodeLayer[i][j] = 0;
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
        NodeLayerPreFru[i] = 0;
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
    for i in 0..20 {
        LeafArea[i] = 0.;
        LeafWeightLayer[i] = 0.;
    }
    LeafWeightLayer[0] = 0.2;
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

fn read_soil_impedance(path: &Path) -> Result<(), Box<dyn std::error::Error>>
//     This function opens the soil root impedance data file and reads it.
//  It is called from ReadInput(), and executed once at the beginning of the
//  simulation. The variables read here are later used to compute soil impedance
//  to root growth.
//
//     Global or file scope variables set:
//        gh2oc, impede, inrim, ncurve, tstbd.
//
{
    let mut rdr = csv::Reader::from_path(path)?;
    unsafe {
        ncurve = 0;
        inrim = 0;
    }
    for result in rdr.records() {
        let record = result?;
        for (i, item) in record.iter().enumerate() {
            if i == 0 {
                let water_content: f64 = item.parse()?;
                unsafe {
                    gh2oc[ncurve as usize] = water_content;
                }
            } else {
                unsafe {
                    impede[i][ncurve as usize] = item.parse()?;
                }
            }
        }
        unsafe {
            ncurve += 1;
        }
    }
    for (j, header) in rdr.headers()?.iter().enumerate() {
        if j > 0 {
            unsafe {
                for i in 0..ncurve as usize {
                    tstbd[j - 1][i] = header.parse()?;
                }
            }
        }
    }
    Ok(())
}
unsafe fn InitSoil(soil_layers: [SoilLayer; 14], soil_hydraulic: SoilHydraulic)
//     This function opens the initial soil data file and reads it. It is
//     executed
//  once at the beginning of the simulation. It is called by ReadInput().
//
//     Global or file scope variables set:
//        rnnh4, rnno3, oma, h2oint.
{
    let mut condfc = [0f64; 9]; // hydraulic conductivity at field capacity of horizon
                                // layers, cm per day.
    let mut h2oint = [0f64; 14]; // initial soil water content, percent of field capacity,
                                 // defined by input for consecutive 15 cm soil layers.
    let mut ldepth = [0f64; 9]; // depth from soil surface to the end of horizon layers, cm.
    let mut oma = [0f64; 14]; // organic matter at the beginning of the season, percent of
                              // soil weight, defined by input for consecutive 15 cm soil
                              // layers.
    let mut pclay = [0f64; 9]; // percentage of clay in soil horizon of horizon layers.
    let mut psand = [0f64; 9]; // percentage of sand in soil horizon of horizon layers.
    let mut psidra: f64; // soil matric water potential, bars, for which immediate
                         // drainage will be simulated (suggested value -0.25 to -0.1).
    let mut psisfc: f64; // soil matric water potential at field capacity,
                         // bars (suggested value -0.33 to -0.1).
    let mut rnnh4 = [0f64; 14]; // residual nitrogen as ammonium in soil at beginning of
                                // season, kg per ha. defined by input for consecutive 15 cm
                                // soil layers.
    let mut rnno3 = [0f64; 14]; // residual nitrogen as nitrate in soil at beginning of
                                // season, kg per ha. defined by input for consecutive 15 cm
                                // soil layers.
    for (i, layer) in soil_layers.iter().enumerate() {
        rnnh4[i] = layer.ammonium;
        rnno3[i] = layer.nitrate;
        oma[i] = layer.organic_matter;
        h2oint[i] = layer.water_content;
    }
    let lyrsol = soil_hydraulic.layers.len();
    //     Zeroise arrays of data values.
    for i in 0..9 {
        airdr[i] = 0.;
        thetas[i] = 0.;
        alpha[i] = 0.;
        beta[i] = 0.;
        SaturatedHydCond[i] = 0.;
        BulkDensity[i] = 0.;
    }
    RatioImplicit = soil_hydraulic.implicit_ratio;
    conmax = soil_hydraulic.max_conductivity;
    psisfc = soil_hydraulic.psi_fc;
    psidra = soil_hydraulic.psi_id;
    if psisfc > 0. {
        psisfc = -psisfc; // make sure it is negative
    }
    if psidra > 0. {
        psidra = -psidra; // make sure it is negative
    }
    //     -Reading next lines  -------------------------------------------
    for (il, layer) in soil_hydraulic.layers.iter().enumerate() {
        //     First line for each layer
        ldepth[il] = layer.depth;
        condfc[il] = layer.hcfc;
        pclay[il] = layer.clay;
        psand[il] = layer.sand;
        airdr[il] = layer.theta_d;
        thetas[il] = layer.theta_s;
        alpha[il] = layer.alpha;
        beta[il] = layer.beta;
        SaturatedHydCond[il] = layer.hcs;
        BulkDensity[il] = layer.bulk_density;
    }
    //
    let mut j = 0usize; // horizon number
    let mut sumdl = 0f64; // depth to the bottom this layer (cm);
    let rm = 2.65f64; // density of the solid fraction of the soil (g / cm3)
    let mut bdl = [0f64; 40]; // array of bulk density of soil layers

    for l in 0..nl as usize {
        //     Using the depth of each horizon layer (ldepth), the horizon
        //  number (SoilHorizonNum) is computed for each soil layer.
        sumdl += dl[l];
        while (sumdl > ldepth[j]) && (j < lyrsol) {
            j += 1;
        }
        SoilHorizonNum[l] = j as i32;
        //     bdl, thad, thts are defined for each soil layer, using the
        //  respective input variables BulkDensity, airdr, thetas.
        //      FieldCapacity, MaxWaterCapacity and thetar are computed for each
        //      layer, as water
        //  content (cm3 cm-3) of each layer corresponding to matric potentials
        //  of psisfc (for field capacity), psidra (for free drainage) and -15
        //  bars (for permanent wilting point), respectively, using function
        //  qpsi.
        //      pore space volume (PoreSpace) is also computed for each layer.
        //      make sure that saturated water content is not more than pore
        //      space.
        bdl[l] = BulkDensity[j];
        PoreSpace[l] = 1. - BulkDensity[j] / rm;
        if thetas[j] > PoreSpace[l] {
            thetas[j] = PoreSpace[l];
        }
        thad[l] = airdr[j];
        thts[l] = thetas[j];
        FieldCapacity[l] = qpsi(psisfc, thad[l], thts[l], alpha[j], beta[j]);
        MaxWaterCapacity[l] = qpsi(psidra, thad[l], thts[l], alpha[j], beta[j]);
        thetar[l] = qpsi(-15., thad[l], thts[l], alpha[j], beta[j]);
        //     When the saturated hydraulic conductivity (SaturatedHydCond) is
        //     not
        //  given, it is computed from the hydraulic conductivity at field
        //  capacity (condfc), using the wcond function.
        if SaturatedHydCond[j] <= 0. {
            SaturatedHydCond[j] =
                condfc[j] / wcond(FieldCapacity[l], thad[l], thts[l], beta[j], 1., 1.);
        }
    }
    //     Loop for all soil layers. Compute depth from soil surface to
    //  the end of each layer (sumdl).
    sumdl = 0.;
    for l in 0..nl as usize {
        sumdl += dl[l];
        //     At start of simulation compute estimated movable fraction of
        //  nitrates in each soil layer, following the work of:
        //     Bowen, W.T., Jones, J.W., Carsky, R.J., and Quintana, J.O. 1993.
        //  Evaluation of the nitrogen submodel of CERES-maize following legume
        //  green manure incorporation. Agron. J. 85:153-159.
        //     The fraction of total nitrate in a layer that is in solution and
        //  can move from one layer to the next with the downward flow of
        //  water, FLOWNO3[l], is a function of the adsorption coefficient,
        //  soil bulk density, and the volumetric soil water content at the
        //  drained upper limit.
        //     Adsorption coefficients are assumed to be 0.0 up to 30 cm depth,
        //  and deeper than 30 cm - 0.2, 0.4, 0.8, 1.0, 1.2, and 1.6 for each
        //  successive 15 cm layer.
        // Adsorption coefficient
        let coeff: f64 = if sumdl <= 30. {
            0.
        } else if sumdl <= 45. {
            0.2
        } else if sumdl <= 60. {
            0.4
        } else if sumdl <= 75. {
            0.6
        } else if sumdl <= 90. {
            0.8
        } else if sumdl <= 105. {
            1.0
        } else if sumdl <= 120. {
            1.2
        } else {
            1.6
        };
        NO3FlowFraction[l] = 1. / (1. + coeff * bdl[l] / MaxWaterCapacity[l]);
        //     Determine the corresponding 15 cm layer of the input file.
        //     Compute the initial volumetric water content (VolWaterContent) of
        //     each
        //  layer, and check that it will not be less than the air-dry value or
        //  more than pore space volume.
        j = ((sumdl - 1.) / 15.).floor() as usize;
        if j > 13 {
            j = 13;
        }
        let n = SoilHorizonNum[l] as usize;
        VolWaterContent[l][0] = FieldCapacity[l] * h2oint[j] / 100.;
        if VolWaterContent[l][0] < airdr[n] {
            VolWaterContent[l][0] = airdr[n];
        }
        if VolWaterContent[l][0] > PoreSpace[l] {
            VolWaterContent[l][0] = PoreSpace[l];
        }
        //     Initial values of ammonium N (rnnh4, VolNh4NContent) and nitrate
        //     N
        //  (rnno3, VolNo3NContent) are converted from kgs per ha to mg / cm3
        //  for each soil layer, after checking for minimal amounts.
        if rnno3[j] < 2.0 {
            rnno3[j] = 2.0;
        }
        if rnnh4[j] < 0.2 {
            rnnh4[j] = 0.2;
        }
        VolNo3NContent[l][0] = rnno3[j] / 15. * 0.01;
        VolNh4NContent[l][0] = rnnh4[j] / 15. * 0.01;
        // organic matter in mg / cm3 units.
        let om = (oma[j] / 100.) * bdl[l] * 1000.;
        //     potom is the proportion of readily mineralizable om. it is a
        //  function of soil depth (sumdl, in cm), modified from GOSSYM (where
        //  it probably includes the 0.4 factor for organic C in om).
        let mut potom = 0.15125 - 0.02878 * sumdl.ln();
        if potom < 0. {
            potom = 0.;
        }
        //     FreshOrganicMatter is the readily mineralizable organic matter (=
        //     "fresh organic
        //  matter" in CERES models). HumusOrganicMatter is the remaining
        //  organic matter, which is mineralized very slowly.
        FreshOrganicMatter[l][0] = om * potom;
        HumusOrganicMatter[l][0] = om * (1. - potom);
    }
    //     Since the initial value has been set for the first column only
    //  in each layer, these values are now assigned to all the other columns.
    for l in 0..nl as usize {
        for k in 1..nk as usize {
            VolWaterContent[l][k] = VolWaterContent[l][0];
            VolNo3NContent[l][k] = VolNo3NContent[l][0];
            VolNh4NContent[l][k] = VolNh4NContent[l][0];
            FreshOrganicMatter[l][k] = FreshOrganicMatter[l][0];
            HumusOrganicMatter[l][k] = HumusOrganicMatter[l][0];
        }
    }
    //     Total amounts of water (InitialTotalSoilWater), nitrate N
    //     (TotalSoilNo3N), ammonium
    //  N (TotalSoilNh4N), and urea N (TotalSoilUreaN) are computed for the
    //  whole slab.
    InitialTotalSoilWater = 0.;
    TotalSoilNo3N = 0.;
    TotalSoilNh4N = 0.;
    TotalSoilUreaN = 0.;
    //
    for l in 0..nl as usize {
        for k in 0..nk as usize {
            InitialTotalSoilWater += VolWaterContent[l][k] * dl[l] * wk[k];
            TotalSoilNo3N += VolNo3NContent[l][k] * dl[l] * wk[k];
            TotalSoilNh4N += VolNh4NContent[l][k] * dl[l] * wk[k];
            VolUreaNContent[l][k] = 0.;
        }
    }
    //     InitialTotalSoilWater is converted from cm3 per slab to mm.
    InitialTotalSoilWater = 10. * InitialTotalSoilWater / RowSpace;
    let bsand = 20f64; // heat conductivity of sand and silt (mcal cm-1 s-1 C-1).
    let bclay = 7f64; // heat conductivity of clay (mcal cm-1 s-1 C-1).
    let cka = 0.0615f64; // heat conductivity of air (mcal cm-1 s-1 C-1).
    let ckw = 1.45f64; // heat conductivity of water (mcal cm-1 s-1 C-1).
    let cmin = 0.46f64; // heat capacity of the mineral fraction of the soil.
    let corg = 0.6f64; // heat capacity of the organic fraction of the soil.
    let ga = 0.144f64; // shape factor for air in pore spaces.
    let ro = 1.3f64; // specific weight of organic fraction of soil.
                     //     Compute aggregation factors:

    dsand = form(bsand, ckw, ga); // aggregation factor for sand in water
    dclay = form(bclay, ckw, ga); // aggregation factor for clay in water
    let dsandair: f64 = form(bsand, cka, ga); // aggregation factor for sand in air
    let dclayair: f64 = form(bclay, cka, ga); // aggregation factor for clay in air
                                              //     Loop over all soil layers, and define indices for some soil arrays.

    sumdl = 0.; // sum of depth of consecutive soil layers.

    for l in 0..nl as usize {
        sumdl += dl[l];
        let mut j = ((sumdl + 14.) / 15.).floor() as usize - 1; //  layer definition for oma
        if j > 13 {
            j = 13;
        }
        //     Using the values of the clay and organic matter percentages in
        //     the soil, compute
        //   mineral and organic fractions of the soil, by weight and by volume.
        let mmo = oma[j] / 100.; // organic matter fraction of dry soil (by weight).
        let mm = 1. - mmo; // mineral fraction of dry soil (by weight).
                           //     MarginalWaterContent is set as a function of the sand fraction of
                           //     the soil.
        let ra = (mmo / ro) / (mm / rm); // volume ratio of organic to mineral soil fractions.

        let i1 = SoilHorizonNum[l] as usize; //  layer definition as in soil hydrology input file.

        //     The volume fractions of clay (ClayVolumeFraction) and of sand
        //     plus silt (SandVolumeFraction), are calculated.
        MarginalWaterContent[l] = 0.1 - 0.07 * psand[i1] / 100.;
        let xo = (1. - PoreSpace[l]) * ra / (1. + ra); // organic fraction of soil (by volume).
        let xm = (1. - PoreSpace[l]) - xo; // mineral fraction of soil (by volume).
        ClayVolumeFraction[l] = pclay[i1] * xm / mm / 100.;
        SandVolumeFraction[l] = 1. - PoreSpace[l] - ClayVolumeFraction[l];
        //     Heat capacity of the solid soil fractions (mineral + organic, by
        //     volume )
        HeatCapacitySoilSolid[l] = xm * cmin + xo * corg;
        //     The heat conductivity of dry soil (HeatCondDrySoil) is computed
        //     using the
        //  procedure suggested by De Vries.
        HeatCondDrySoil[l] = 1.25
            * (PoreSpace[l] * cka
                + dsandair * bsand * SandVolumeFraction[l]
                + dclayair * bclay * ClayVolumeFraction[l])
            / (PoreSpace[l] + dsandair * SandVolumeFraction[l] + dclayair * ClayVolumeFraction[l]);
    }
}

unsafe fn ReadPlantMapInput(plant_maps: Vec<PlantMap>)
//     This sunbroutine opens and reads an ascii file with input of
//  observed plant map adjustment data. It is used to adjust the simulation.
//     It is called by ReadInput().
//     The following global variables are referenced:    iyear, PlantmapFileName
//     The following global variables are set:
//       MapDataGreenBollNum[], MapDataDate[], MapDataMainStemNodes[],
//       MapDataPlantHeight[], MapDataSquareNum[], MapDataAllSiteNum[];
//
{
    for (i, plant_map) in plant_maps.iter().enumerate() {
        MapDataDate[i] = plant_map.date.ordinal() as i32; // day of year
        MapDataPlantHeight[i] = plant_map.plant_height; // Plant height, cm
        MapDataMainStemNodes[i] = plant_map.main_stem_nodes; // Number of mainstem nodes
        MapDataSquareNum[i] = plant_map.number_of_squares; // Number of squares per plant
        MapDataGreenBollNum[i] = plant_map.number_of_bolls; // Number of green bolls per plant
        MapDataAllSiteNum[i] = plant_map.number_of_nodes; // Number of total sites per plant
    }
}
