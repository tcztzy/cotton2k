use crate::meteorology::Meteorology;
use crate::plant_growth::{PhysiologicalAge, PlantGrowth};
use crate::plant_nitrogen::PlantNitrogen;
use crate::soil_temperature::SoilThermology;
use crate::utils::{cell_distance, fmax, fmin, slab_horizontal_location, slab_vertical_location};
use crate::{
    addwtbl, beta, dl, isw, light_intercept_parameters, maxl, nk, nl, noitr, pixday, rracol, thad,
    thts, wcond, wk, ActualTranspiration, AgeOfPreFruNode, AgronomyOperation, AppliedWater,
    AverageLeafAge, AverageLwp, AverageLwpMin, AveragePsi, AverageSoilPsi, BurrWeightOpenBolls,
    CapillaryFlow, CheckDryMatterBal, ComputeIrrigation, Cotton2KError, CottonPhenology,
    CottonWeightOpenBolls, CumFertilizerN, CumNetPhotosynth, CumNitrogenUptake, CumTranspiration,
    CumWaterAdded, CumWaterDrained, DayEmerge, DayInc, DayLength, DayOfSimulation, DayStart,
    DayStartPredIrrig, DayStopPredIrrig, DayTimeTemp, Daynum, Defoliate, Drain, ElCondSatSoilToday,
    FertilizationMethod, FirstSquare, GetFromClim, Irrig, IrrigMethod, Kday, LeafAge, LeafArea,
    LeafAreaIndex, LeafAreaIndexes, LeafAreaMainStem, LeafAreaNodes, LeafAreaPreFru, LeafNConc,
    LeafNitrogen, LeafResistance, LightIntercept, LightInterceptLayer, LightInterceptMethod,
    LocationColumnDrip, LocationLayerDrip, LwpMax, LwpMin, LwpMinX, LwpX, MaxIrrigation,
    MaxWaterCapacity, NO3FlowFraction, NetPhotosynthesis, NodeLayer, NodeLayerPreFru,
    NumFruitBranches, NumIrrigations, NumLayersWithRoots, NumNodes, NumPreFruNodes, NumVegBranches,
    PerPlantArea, PlantHeight, PlantPopulation, PlantRowColumn, PlantWeight, PoreSpace, Profile,
    ReferenceETP, RootColNumLeft, RootWtCapblUptake, RootsCapableOfUptake, RowSpace,
    SaturatedHydCond, Scratch21, SoilNitrogen, SoilNitrogenAverage, SoilNitrogenBal,
    SoilNitrogenLoss, SoilPsi, SoilSum, StemWeight, SupplyNH4N, SupplyNO3N, TotalLeafWeight,
    VolNh4NContent, VolNo3NContent, VolUreaNContent, VolWaterContent, WaterStress, WaterStressStem,
    WaterTableLayer, WaterUptake, CLIMATE_METRIC_IRRD, CLIMATE_METRIC_RAIN,
};
use chrono::{Datelike, NaiveDate};

#[derive(Debug, Clone, Copy)]
pub struct State {
    pub date: NaiveDate,

    pub plant_height: f64,
    pub rracol: [f64; 20],
    /// residual available carbon for root growth from previous day.
    pub pavail: f64,

    // for plant_nitrogen
    /// daily added nitrogen to fruit, g per plant.
    pub addnf: f64,
    /// daily added nitrogen to root, g per plant.
    pub addnr: f64,
    /// daily added nitrogen to vegetative shoot, g per plant.
    pub addnv: f64,
    /// amount of nitrogen not used for growth of plant parts.
    pub xtran: f64,
    /// reserve N in leaves, in g per plant.
    pub leafrs: f64,
    /// reserve N in petioles, in g per plant.
    pub petrs: f64,
    /// reserve N in stems, in g per plant.
    pub stemrs: f64,
    /// reserve N in roots, in g per plant.
    pub rootrs: f64,
    /// reserve N in burrs, in g per plant.
    pub burres: f64,
    /// nitrogen requirement for fruit growth.
    pub reqf: f64,
    /// total nitrogen requirement for plant growth.
    pub reqtot: f64,
    /// nitrogen requirement for vegetative shoot growth.
    pub reqv: f64,
    /// nitrogen requirement for burr growth.
    pub rqnbur: f64,
    /// nitrogen requirement for leaf growth.
    pub rqnlef: f64,
    /// nitrogen requirement for petiole growth.
    pub rqnpet: f64,
    /// nitrogen requirement for root growth.
    pub rqnrut: f64,
    /// nitrogen requirement for seed growth.
    pub rqnsed: f64,
    /// nitrogen requirement for square growth.
    pub rqnsqr: f64,
    /// nitrogen requirement for stem growth.
    pub rqnstm: f64,
    /// total nitrogen available for growth.
    pub npool: f64,
    /// nitrogen uptake from the soil, g per plant.
    pub uptn: f64,
}

impl State {
    pub fn new(date: NaiveDate) -> Self {
        State {
            date,
            plant_height: 4.0,
            rracol: [1.; 20],
            pavail: 0.,
            addnf: 0.,
            addnr: 0.,
            addnv: 0.,
            xtran: 0.,
            leafrs: 0.,
            petrs: 0.,
            stemrs: 0.,
            rootrs: 0.,
            burres: 0.,
            reqf: 0.,
            reqtot: 0.,
            reqv: 0.,
            rqnbur: 0.,
            rqnlef: 0.,
            rqnpet: 0.,
            rqnrut: 0.,
            rqnsed: 0.,
            rqnsqr: 0.,
            rqnstm: 0.,
            npool: 0.,
            uptn: 0.,
        }
    }

    /// This function executes all the simulation computations in a day. It is called from [Profile::run()], and
    /// [Profile::adjust()].
    ///
    /// It calls the following functions:
    /// * [Profile::column_shading()]
    /// * [soil_thermology()]
    /// * [SoilProcedures()]
    /// * [SoilNitrogen()]
    /// * [SoilSum()]
    /// * [PhysiologicalAge()]
    /// * [Profile::pix()]
    /// * [Defoliate()]
    /// * [Profile::stress()]
    /// * [Profile::get_net_photosynthesis()]
    /// * [PlantGrowth::plant_growth()]
    /// * [CottonPhenology()]
    /// * [PlantNitrogen::plant_nitrogen()]
    /// * [CheckDryMatterBal()]
    /// * [PlantNitrogen::plant_nitrogen_balance()]
    /// * [SoilNitrogenBal()]
    /// * [SoilNitrogenAverage()]
    ///
    /// The following global variables are referenced here:
    /// * [DayEmerge]
    /// * [DayFinish]
    /// * [DayStart]
    /// * [Kday]
    /// * [LeafAreaIndex]
    /// * [pixday]
    ///
    /// The following global variables are set here:
    /// * [bEnd]
    /// * [DayInc]
    /// * [Daynum]
    /// * [DayOfSimulation]
    /// * [isw]
    /// * [Kday]
    pub fn simulate_this_day(&mut self, profile: &mut Profile) -> Result<(), Cotton2KError> {
        unsafe {
            // Compute Daynum (day of year), Date, and DayOfSimulation (days from start of simulation).
            Daynum += 1;
            DayOfSimulation = Daynum - DayStart + 1;
            // Compute Kday (days from emergence).
            Kday = if DayEmerge <= 0 {
                0
            } else {
                Daynum - DayEmerge + 1
            };
            if Kday < 0 {
                Kday = 0;
            }
            // The following functions are executed each day (also before emergence).
            self.column_shading(profile); // computes light interception and soil shading.
            self.meteorology(profile); // computes climate variables for today.
            self.soil_thermology(); // executes all modules of soil and canopy temperature.
            self.soil_procedures(profile)?; // executes all other soil processes.
            SoilNitrogen(); // computes nitrogen transformations in the soil.
            SoilSum(); // computes totals of water and N in the soil.

            // The following is executed each day after plant emergence:
            if Daynum >= DayEmerge && isw > 0 {
                // If this day is after emergence, assign to isw the value of 2.
                isw = 2;
                DayInc = PhysiologicalAge(); // computes physiological age
                if pixday[0] > 0 {
                    self.pix(); // effects of pix applied.
                }
                Defoliate(); // effects of defoliants applied.
                self.stress(profile); // computes water stress factors.
                self.get_net_photosynthesis(profile)?; // computes net photosynthesis.
                self.plant_growth(); // executes all modules of plant growth.
                CottonPhenology(); // executes all modules of plant phenology.
                self.plant_nitrogen(); // computes plant nitrogen allocation.
                CheckDryMatterBal(); // checks plant dry matter balance.

                // If the relevant output flag is not zero, compute soil nitrogen balance and soil nitrogen averages by
                // layer, and write this information to files.
                if false {
                    self.plant_nitrogen_balance(); // checks plant nitrogen balance.
                    SoilNitrogenBal(); // checks soil nitrogen balance.
                    SoilNitrogenAverage(); // computes average soil nitrogen by layers.
                }
            }
            // Check if the date to stop simulation has been reached, or if this is the last day with available weather
            // data. Simulation will also stop when no leaves remain on the plant.
            if self.date >= profile.last_day_weather_data {
                return Err(Cotton2KError {
                    level: 0,
                    message: String::from("No more weather data!"),
                });
            }
            if Kday > 10 && LeafAreaIndex < 0.0002 {
                return Err(Cotton2KError {
                    level: 0,
                    message: String::from("Leaf area index is too small!"),
                });
            }
        }
        Ok(())
    }
    /// This function computes light interception by crop canopy and shading of soil columns by the plants. It is called from SimulateThisDay().
    ///
    /// The following global variables are referenced here:
    /// * [DayEmerge]
    /// * [Daynum]
    /// * [isw]
    /// * [LeafAreaIndex]
    /// * [PlantHeight]
    /// * [PlantRowColumn]
    /// * [nk]
    /// * [RowSpace]
    ///
    /// The following global variables are set here:
    /// * [LightIntercept]
    /// * [rracol]
    fn column_shading(&mut self, profile: &mut Profile) {
        unsafe {
            // Before emergence: no light interception and no shading.
            // LightIntercept is assigned zero, and the rracol array is assigned 1.
            if Daynum < DayEmerge || isw <= 0 || DayEmerge <= 0 {
                LightIntercept = 0.;
                for k in 0..nk as usize {
                    rracol[k] = 1.;
                }
                return;
            }
            // Compute the maximum leaf area index until this day (lmax).
            if LeafAreaIndex > profile.lmax {
                profile.lmax = LeafAreaIndex;
            }
            // Light interception is computed by two methods:
            //
            // 1. It is assumed to be proportional to the ratio of plant height to row spacing.
        }
        // light interception computed from plant height.
        let zint = unsafe { 1.0756 * PlantHeight / RowSpace };
        match profile.light_intercept_method {
            LightInterceptMethod::Latered => unsafe {
                for i in 0..20 {
                    LeafArea[i] = 0.;
                    AverageLeafAge[i] = 0.;
                }
                for i in 0..9 {
                    LeafArea[NodeLayerPreFru[i] as usize] += LeafAreaPreFru[i];
                    AverageLeafAge[NodeLayerPreFru[i] as usize] +=
                        LeafAreaPreFru[i] * AgeOfPreFruNode[i];
                }
                for k in 0..NumVegBranches as usize {
                    for l in 0..NumFruitBranches[k] as usize {
                        LeafArea[NodeLayer[k][l] as usize] += LeafAreaMainStem[k][l];
                        AverageLeafAge[NodeLayer[k][l] as usize] +=
                            LeafAreaMainStem[k][l] * LeafAge[k][l][0];
                        for m in 0..NumNodes[k][l] as usize {
                            LeafArea[NodeLayer[k][l] as usize] += LeafAreaNodes[k][l][m];
                            AverageLeafAge[NodeLayer[k][l] as usize] +=
                                LeafAreaNodes[k][l][m] * LeafAge[k][l][m];
                        }
                    }
                }
                if FirstSquare <= 0 {
                    LeafArea[0] += 0.20 * 0.6;
                }
                let mut light_through = 0.;
                for i in 0..20 {
                    AverageLeafAge[i] /= LeafArea[i];
                    LeafAreaIndexes[i] = LeafArea[i] / PerPlantArea;
                    LightInterceptLayer[i] =
                        1. - (light_intercept_parameters[i] * LeafAreaIndexes[i]).exp();
                    light_through += light_intercept_parameters[i] * LeafAreaIndexes[i];
                }
                LightIntercept = 1. - light_through.exp();
            },
            LightInterceptMethod::Fry1980 => unsafe {
                LightIntercept = 0.39 * LeafAreaIndex.powf(0.68);
            },
            LightInterceptMethod::Original => {
                unsafe {
                    // 2. It is computed as a function of leaf area index. If LeafAreaIndex is not greater than 0.5 lfint is a linear function of it.

                    // light interception computed from leaf area index.
                    let lfint = if LeafAreaIndex <= 0.5 {
                        0.80 * LeafAreaIndex
                    } else {
                        // If the leaf area index is greater than 0.5, lfint is computed as an exponential function of LeafAreaIndex.
                        1. - (0.07 - 1.16 * LeafAreaIndex).exp()
                    };
                    // If lfint is greater then zint, LightIntercept is their average value.
                    // Otherwise, if the LeafAreaIndex is decreasing, it is lfint. Else it is zint.
                    LightIntercept = if lfint > zint {
                        0.5 * (zint + lfint)
                    } else if LeafAreaIndex < profile.lmax {
                        lfint
                    } else {
                        zint
                    };
                }
            }
        }
        unsafe {
            // The value of LightIntercept is between zero and one.
            if LightIntercept < 0. {
                LightIntercept = 0.;
            }
            if LightIntercept > 1. {
                LightIntercept = 1.;
            }
        }
        // Loop of soil columns.
        let mut sw = 0.; // sum of column widths
        let mut sw0: f64 = 0.; // sum of column widths up to location of plant row.
        let mut sw1: f64; // distance from middle of a column to the plant row, cm.
        let mut j;
        let mut k0; // number of columns from plant row location.
        unsafe {
            for k in 0..nk as usize {
                if k <= PlantRowColumn as usize {
                    // When the column is on the left of the plant row.
                    j = (PlantRowColumn - k as i32) as usize;
                    sw += wk[j];
                    sw0 = sw;
                    sw1 = sw - wk[j] / 2.;
                    k0 = j;
                } else {
                    // When the column is on the right of the plant row.
                    sw += wk[k];
                    sw1 = sw - sw0 - wk[k] / 2.;
                    k0 = k;
                }

                // Relative shading is computed as a function of sw1 and plant height, modified by LightIntercept.
                // [rracol] is the fraction of radiation received by a soil column. Iys minimum value is 0.05 .

                // relative shading of a soil column.
                let shade = if sw1 >= PlantHeight {
                    0.
                } else {
                    let mut result = 1. - (sw1 / PlantHeight).powi(2);
                    if LightIntercept < zint && LeafAreaIndex < profile.lmax {
                        result *= LightIntercept / zint;
                    }
                    result
                };
                rracol[k0] = 1. - shade;
                if rracol[k0] < 0.05 {
                    rracol[k0] = 0.05;
                }
            }
        }
    }
    /// This function manages all the soil related processes, and is executed once each day.
    ///
    /// It is called from [Profile::simulate_this_day()] and it calls the following functions:
    /// * [ApplyFertilizer()]
    /// * [AveragePsi()]
    /// * [CapillaryFlow()]
    /// * [ComputeIrrigation()]
    /// * [DripFlow()]
    /// * [GravityFlow()]
    /// * [RootsCapableOfUptake()]
    /// * [WaterUptake()]
    /// * [WaterTable()]
    ///
    /// The following global variables are referenced here:
    /// * [ActualTranspiration]
    /// * [Clim]
    /// * [DayEmerge]
    /// * [Daynum]
    /// * [DayStartPredIrrig]
    /// * [DayStopPredIrrig]
    /// * [dl]
    /// * [Irrig]
    /// * [IrrigMethod]
    /// * [isw]
    /// * [Kday]
    /// * [MaxIrrigation]
    /// * [nk]
    /// * [nl]
    /// * [NumIrrigations]
    /// * [PerPlantArea]
    /// * [SoilPsi]
    /// * [RowSpace]
    /// * [SupplyNH4N]
    /// * [SupplyNO3N]
    /// * [VolWaterContent]
    /// * [wk]
    ///
    /// The following global variables are set here:
    /// * [AverageSoilPsi]
    /// * [CumNitrogenUptake]
    /// * [CumTranspiration]
    /// * [CumWaterAdded]
    /// * [LocationColumnDrip]
    /// * [LocationLayerDrip]
    /// * [noitr]
    fn soil_procedures(&mut self, profile: &mut Profile) -> Result<(), Cotton2KError> {
        // The following constant parameters are used:
        const cpardrip: f64 = 0.2;
        const cparelse: f64 = 0.4;
        // Call function ApplyFertilizer() for nitrogen fertilizer application.
        self.apply_fertilizer(profile)?;
        let mut drip_water_amount = 0f64; // amount of water applied by drip irrigation

        // amount of water applied by non-drip irrigation or rainfall
        // Check if there is rain on this day
        let mut water_to_apply = unsafe { GetFromClim(CLIMATE_METRIC_RAIN, Daynum) };
        // If irrigation is to be predicted for this day, call ComputeIrrigation() to compute the actual amount of irrigation.
        unsafe {
            if MaxIrrigation > 0. {
                if Daynum >= DayStartPredIrrig && Daynum < DayStopPredIrrig {
                    ComputeIrrigation();
                    if IrrigMethod == 2 {
                        drip_water_amount = AppliedWater;
                    } else {
                        water_to_apply += AppliedWater;
                    }
                    AppliedWater = 0.;
                }
            }
        }
        // When water is added by an irrigation defined in the input: update the amount of applied water.
        unsafe {
            for i in 0..NumIrrigations as usize {
                if Daynum == Irrig[i].day {
                    if Irrig[i].method == 2 {
                        drip_water_amount += Irrig[i].amount;
                        LocationColumnDrip = Irrig[i].LocationColumnDrip;
                        LocationLayerDrip = Irrig[i].LocationLayerDrip;
                    } else {
                        water_to_apply += Irrig[i].amount;
                    }
                    break;
                }
            }
        }
        unsafe {
            CumWaterAdded += water_to_apply + drip_water_amount;
        }
        unsafe {
            if Kday > 0 {
                Scratch21[(DayOfSimulation - 1) as usize].amitri =
                    water_to_apply + drip_water_amount;
            }
        }
        unsafe {
            // The following will be executed only after plant emergence
            if Daynum >= DayEmerge && isw > 0 {
                RootsCapableOfUptake(); // function computes roots capable of uptake for each soil cell
                AverageSoilPsi = AveragePsi(); // function computes the average matric soil water potential in the root zone, weighted by the roots-capable-of-uptake.
                WaterUptake(); // function  computes water and nitrogen uptake by plants.

                // Update the cumulative sums of actual transpiration (CumTranspiration, mm) and total uptake of nitrogen (CumNitrogenUptake, mg N per slab, converted from total N supply, g per plant).
                CumTranspiration += ActualTranspiration;
                CumNitrogenUptake += (SupplyNO3N + SupplyNH4N) * 10. * RowSpace / PerPlantArea;
            }
        }
        // Call function WaterTable() for saturating soil below water table.
        if profile.num_watertable_data > 0 {
            self.watertable(profile)?;
        }
        if water_to_apply > 0. {
            // For rain or surface irrigation.
            // The number of iterations is computed from the thickness of the first soil layer.
            unsafe {
                noitr = (cparelse * water_to_apply / (dl[0] + 2.) + 1.) as i32;
            }
            // the amount of water applied, mm per iteration.
            let applywat = water_to_apply / unsafe { noitr } as f64;
            // The following subroutines are called noitr times per day:
            // If water is applied, GravityFlow() is called when the method of irrigation is not by drippers, followed by CapillaryFlow().
            for _ in 0..unsafe { noitr } as usize {
                unsafe {
                    GravityFlow(applywat);
                    CapillaryFlow();
                }
            }
        }
        if drip_water_amount > 0. {
            // For drip irrigation.
            // The number of iterations is computed from the volume of the soil cell in which the water is applied.
            unsafe {
                noitr = (cpardrip * drip_water_amount
                    / (dl[LocationLayerDrip as usize] * wk[LocationColumnDrip as usize])
                    + 1.) as i32;
            }
            // the amount of water applied, mm per iteration.
            let applywat = drip_water_amount / unsafe { noitr } as f64;
            // If water is applied, DripFlow() is called followed by CapillaryFlow().
            for _ in 0..unsafe { noitr } as usize {
                unsafe {
                    DripFlow(applywat)?;
                    CapillaryFlow();
                }
            }
        }
        // When no water is added, there is only one iteration in this day.
        if water_to_apply + drip_water_amount <= 0. {
            unsafe {
                noitr = 1;
                CapillaryFlow();
            }
        }
        Ok(())
    }

    /// This function simulates the application of nitrogen fertilizer on each date of application.
    /// It is called from [SoilProcedures()].
    ///
    /// The following global variables are referenced here:
    /// * [Daynum]
    /// * [dl]
    /// * [LightIntercept]
    /// * [LocationColumnDrip]
    /// * [LocationLayerDrip]
    /// * [nk]
    /// * [nl]
    /// * [RowSpace]
    /// * [wk]
    ///
    /// The following global variables are set here:
    /// * [CumFertilizerN]
    /// * [LeafNitrogen]
    /// * [VolNh4NContent]
    /// * [VolNo3NContent]
    /// * [VolUreaNContent]
    fn apply_fertilizer(&mut self, profile: &Profile) -> Result<(), Cotton2KError> {
        let ferc = 0.01; // constant used to convert kgs per ha to mg cm-2

        for ao in &profile.agronomy_operations {
            match ao {
                AgronomyOperation::fertilization {
                    date,
                    urea,
                    nitrate,
                    ammonium,
                    method,
                    drip_x,
                    drip_y,
                } => {
                    unsafe {
                        if Daynum as u32 == date.ordinal() {
                            CumFertilizerN += ferc * RowSpace * (ammonium + nitrate + urea);
                            match method {
                                FertilizationMethod::Broadcast => {
                                    // Compute the number of layers affected by broadcast fertilizer incorporation (lplow), assuming that the depth of incorporation is 20 cm.
                                    let mut lplow: usize = 0; // number of soil layers affected by cultivation
                                    let mut sdl: f64 = 0.; // sum of depth of consecutive soil layers
                                    for l in 0..nl as usize {
                                        sdl += dl[l];
                                        if sdl >= 20. {
                                            lplow = l + 1;
                                            break;
                                        }
                                    }
                                    // Calculate the actual depth of fertilizer incorporation in the soil (fertdp) as the sum of all soil layers affected by incorporation.
                                    let mut fertdp = 0f64; // depth of broadcast fertilizer incorporation, cm
                                    for l in 0..lplow {
                                        fertdp += dl[l];
                                    }
                                    // Update the nitrogen contents of all soil soil cells affected by this fertilizer application.
                                    for l in 0..lplow {
                                        for k in 0..nk as usize {
                                            VolNh4NContent[l][k] += ammonium * ferc / fertdp;
                                            VolNo3NContent[l][k] += nitrate * ferc / fertdp;
                                            VolUreaNContent[l][k] += urea * ferc / fertdp;
                                        }
                                    }
                                }
                                FertilizationMethod::Foliar => {
                                    // It is assumed that 70% of the amount of ammonium or urea intercepted by the canopy is added to the leaf N content (LeafNitrogen).
                                    LeafNitrogen +=
                                        0.70 * LightIntercept * (ammonium + urea) * 1000.
                                            / PlantPopulation;
                                    // The amount not intercepted by the canopy is added to the soil. If the fertilizer is nitrate, it is assumed that all of it is added to the upper soil layer.
                                    // Update nitrogen contents of the upper layer.
                                    for k in 0..nk as usize {
                                        VolNh4NContent[0][k] +=
                                            ammonium * (1. - 0.70 * LightIntercept) * ferc / dl[0];
                                        VolNo3NContent[0][k] += nitrate * ferc / dl[0];
                                        VolUreaNContent[0][k] +=
                                            urea * (1. - 0.70 * LightIntercept) * ferc / dl[0];
                                    }
                                }
                                FertilizationMethod::Sidedress => {
                                    // Define the soil column (ksdr) and the soil layer (lsdr) in which the side-dressed fertilizer is applied.
                                    let ksdr = slab_horizontal_location(*drip_x, RowSpace)?; // the column in which the side-dressed is applied
                                    let lsdr = slab_vertical_location(*drip_y)?; // the layer in which the side-dressed is applied
                                    let mut n00: u32 = 1; // number of soil soil cells in which side-dressed fertilizer is incorporated.

                                    // If the volume of this soil cell is less than 100 cm3, it is assumed that the fertilizer is also incorporated in the soil cells below and to the sides of it.
                                    if (dl[lsdr] * wk[ksdr]) < 100. {
                                        if ksdr < (nk - 1) as usize {
                                            n00 += 1;
                                        }
                                        if ksdr > 0 {
                                            n00 += 1;
                                        }
                                        if lsdr < (nl - 1) as usize {
                                            n00 += 1;
                                        }
                                    }
                                    // amount of ammonium N added to the soil by sidedressing (mg per cell)
                                    let addamm = ammonium * ferc * RowSpace / n00 as f64;
                                    // amount of nitrate N added to the soil by sidedressing (mg per cell)
                                    let addnit = nitrate * ferc * RowSpace / n00 as f64;
                                    // amount of urea N added to the soil by  sidedressing (mg per cell)
                                    let addnur = urea * ferc * RowSpace / n00 as f64;
                                    // Update the nitrogen contents of these soil cells.
                                    VolNo3NContent[lsdr][ksdr] += addnit / (dl[lsdr] * wk[ksdr]);
                                    VolNh4NContent[lsdr][ksdr] += addamm / (dl[lsdr] * wk[ksdr]);
                                    VolUreaNContent[lsdr][ksdr] += addnur / (dl[lsdr] * wk[ksdr]);
                                    if (dl[lsdr] * wk[ksdr]) < 100. {
                                        if ksdr < (nk - 1) as usize {
                                            let kp1 = ksdr + 1; // column to the right of ksdr.
                                            VolNo3NContent[lsdr][kp1] +=
                                                addnit / (dl[lsdr] * wk[kp1]);
                                            VolNh4NContent[lsdr][kp1] +=
                                                addamm / (dl[lsdr] * wk[kp1]);
                                            VolUreaNContent[lsdr][kp1] +=
                                                addnur / (dl[lsdr] * wk[kp1]);
                                        }
                                        if ksdr > 0 {
                                            let km1 = ksdr - 1; // column to the left of ksdr.
                                            VolNo3NContent[lsdr][km1] +=
                                                addnit / (dl[lsdr] * wk[km1]);
                                            VolNh4NContent[lsdr][km1] +=
                                                addamm / (dl[lsdr] * wk[km1]);
                                            VolUreaNContent[lsdr][km1] +=
                                                addnur / (dl[lsdr] * wk[km1]);
                                        }
                                        if lsdr < (nl - 1) as usize {
                                            let lp1 = lsdr + 1;
                                            VolNo3NContent[lp1][ksdr] +=
                                                addnit / (dl[lp1] * wk[ksdr]);
                                            VolNh4NContent[lp1][ksdr] +=
                                                addamm / (dl[lp1] * wk[ksdr]);
                                            VolUreaNContent[lp1][ksdr] +=
                                                addnur / (dl[lp1] * wk[ksdr]);
                                        }
                                    }
                                }
                                FertilizationMethod::Drip => {
                                    // Convert amounts added to mg cm-3, and update the nitrogen content of the soil cell in which the drip outlet is situated.
                                    VolNh4NContent[LocationLayerDrip as usize]
                                        [LocationColumnDrip as usize] += ammonium * ferc * RowSpace
                                        / (dl[LocationLayerDrip as usize]
                                            * wk[LocationColumnDrip as usize]);
                                    VolNo3NContent[LocationLayerDrip as usize]
                                        [LocationColumnDrip as usize] += nitrate * ferc * RowSpace
                                        / (dl[LocationLayerDrip as usize]
                                            * wk[LocationColumnDrip as usize]);
                                    VolUreaNContent[LocationLayerDrip as usize]
                                        [LocationColumnDrip as usize] += urea * ferc * RowSpace
                                        / (dl[LocationLayerDrip as usize]
                                            * wk[LocationColumnDrip as usize]);
                                }
                            }
                        }
                    }
                }
                _ => {}
            }
        }
        Ok(())
    }

    /// effects of pix applied.
    ///
    /// TODO
    fn pix(self: &Self) {}

    /// This function computes the water stress variables affecting the cotton plants.
    /// It is called by [Profile::simulate_this_day()] and calls [LeafWaterPotential()].
    ///
    /// The following global variables are referenced here:
    /// * [Kday]
    /// * [LwpMin]
    /// * [LwpMax]
    /// The following global variables are set here:
    /// * [AverageLwp]
    /// * [AverageLwpMin]
    /// * [LwpMinX]
    /// * [LwpX]
    /// * [Profile::ptsred]
    /// * [WaterStress]
    /// * [WaterStressStem]
    unsafe fn stress(&mut self, profile: &mut Profile) {
        // The following constant parameters are used:
        const vstrs: [f64; 9] = [-3.0, 3.229, 1.907, 0.321, -0.10, 1.230, 0.340, 0.30, 0.05];
        // Call LeafWaterPotential() to compute leaf water potentials.
        LeafWaterPotential();
        // The running averages, for the last three days, are computed: AverageLwpMin is the average of LwpMin, and
        // AverageLwp of LwpMin + LwpMax.
        AverageLwpMin += (LwpMin - LwpMinX[2]) / 3.;
        AverageLwp += (LwpMin + LwpMax - LwpX[2]) / 3.;
        for i in [2, 1] {
            LwpMinX[i] = LwpMinX[i - 1];
            LwpX[i] = LwpX[i - 1];
        }
        LwpMinX[0] = LwpMin;
        LwpX[0] = LwpMin + LwpMax;
        //     No stress effects before 5th day after emergence.
        if Kday < 5 {
            profile.ptsred = 1.;
            WaterStress = 1.;
            WaterStressStem = 1.;
            return;
        }
        // The computation of ptsred, the effect of moisture stress on the photosynthetic rate, is based on the
        // following work:
        // * Ephrath, J.E., Marani, A., Bravdo, B.A., 1990. Effects of moisture stress on stomatal resistance and
        //   photosynthetic rate in cotton (Gossypium hirsutum) 1. Controlled levels of stress. Field Crops Res.
        //   23:117-131.
        //
        // It is a function of AverageLwpMin (average LwpMin for the last three days).
        if AverageLwpMin < vstrs[0] {
            AverageLwpMin = vstrs[0];
        }
        profile.ptsred = vstrs[1] + AverageLwpMin * (vstrs[2] + vstrs[3] * AverageLwpMin);
        if profile.ptsred > 1. {
            profile.ptsred = 1.;
        }
        // The general moisture stress factor (WaterStress) is computed as an empirical function of AverageLwp. psilim,
        // the value of AverageLwp at the maximum value of the function, is used for truncating it.

        // The minimum value of WaterStress is 0.05, and the maximum is 1.
        let psilim = -0.5 * vstrs[5] / vstrs[6]; // limiting value of AverageLwp.
        if AverageLwp > psilim {
            WaterStress = 1.;
        } else {
            WaterStress = vstrs[4] - AverageLwp * (vstrs[5] + vstrs[6] * AverageLwp);
            if WaterStress > 1. {
                WaterStress = 1.;
            }
            if WaterStress < 0.05 {
                WaterStress = 0.05;
            }
        }
        // Water stress affecting plant height and stem growth (WaterStressStem) is assumed to be more severe than
        // WaterStress, especially at low WaterStress values.
        WaterStressStem = WaterStress * (1. + vstrs[7] * (2. - WaterStress)) - vstrs[7];
        if WaterStressStem < vstrs[8] {
            WaterStressStem = vstrs[8];
        }
    }
    /// This function simulates the net photosynthesis of cotton  plants. It is called daily by
    /// [Profile::simulate_this_day()]. This is essentially the routine of GOSSYM with minor changes.
    ///
    /// The following global and file scope variables are referenced here:
    /// * [BurrWeightOpenBolls]
    /// * [CottonWeightOpenBolls]
    /// * [DayLength]
    /// * [Daynum]
    /// * [DayEmerge]
    /// * [DayTimeTemp]
    /// * [iyear]
    /// * [Kday]
    /// * [LeafNConc]
    /// * [LightIntercept]
    /// * [PerPlantArea]
    /// * [PlantWeight]
    /// * [Profile::ptsred]
    /// * [StemWeight]
    /// * [TotalLeafWeight()]
    ///
    /// The following global variables are set here:
    /// * [bEnd]
    /// * [CumNetPhotosynth]
    /// * [NetPhotosynthesis]
    ///
    /// References:
    ///
    /// * Baker et. al. (1972). Simulation of Growth and Yield in Cotton: I. Gross photosynthesis, respiration and
    ///   growth. Crop Sci. 12:431-435.
    /// * Harper et. al. (1973) Carbon dioxide and the photosynthesis of field crops. A metered carbon dioxide release
    ///   in cotton under field conditions.  Agron. J. 65:7-11.
    /// * Baker (1965) Effects of certain environmental factors on net assimilation in cotton. Crop Sci. 5:53-56 (Fig 5).
    unsafe fn get_net_photosynthesis(&mut self, profile: &Profile) -> Result<(), Cotton2KError> {
        //  constants:
        const gsubr: f64 = 0.375; // the growth resiration factor.
        const rsubo: f64 = 0.0032; // maintenance respiration factor.
        const vpnet: [f64; 4] = [1.30, 0.034, 0.010, 0.32];
        const co2parm: [f64; 45] =
            // parameters used to correct photosynthesis for ambient CO2 concentration.
            [
                1.0235, 1.0264, 1.0285, 1.0321, 1.0335, 1.0353, 1.0385, 1.0403, 1.0431, 1.0485,
                1.0538, 1.0595, 1.0627, 1.0663, 1.0716, 1.0752, 1.0784, 1.0823, 1.0880, 1.0923,
                1.0968, 1.1019, 1.1087, 1.1172, 1.1208, 1.1243, 1.1311, 1.1379, 1.1435, 1.1490,
                1.1545, 1.1601, 1.1656, 1.1712, 1.1767, 1.1823, 1.1878, 1.1934, 1.1990, 1.2045,
                1.2101, 1.2156, 1.2212, 1.2267, 1.2323,
            ];
        // Note: co2parm is for icrease in ambient CO2 concentration changes from 1959 (308 ppm).
        // The first 28 values (up to 1987) are from GOSSYM. The other values (up to 2004) are derived from data of the
        // Carbon Dioxide Information Analysis Center (CDIAC).
        //
        // Exit the function and end simulation if there are no leaves.

        if TotalLeafWeight() <= 0. {
            return Err(Cotton2KError {
                level: 0,
                message: String::from("Leaf weight is less than 0!"),
            });
        }

        // Get the CO2 correction factor (pnetcor) for photosynthesis, using AmbientCO2Factor and a factor that may be
        // variety specific (vpnet[0]).

        // correction factor for gross photosynthesis.
        let pnetcor = profile.ambient_CO2_factor
            * vpnet[0]
            * profile.co2_enrichment.unwrap_or_default().factor;
        // Compute ptnfac, the effect of leaf N concentration on photosynthesis, using an empirical relationship.
        // correction factor for low nitrogen content in leaves.
        let mut ptnfac =
            vpnet[3] + (LeafNConc - vpnet[2]) * (1. - vpnet[3]) / (vpnet[1] - vpnet[2]);
        if ptnfac > 1. {
            ptnfac = 1.;
        }
        if ptnfac < vpnet[3] {
            ptnfac = vpnet[3];
        }
        // Convert the average daily short wave radiation from langley per day, to Watts per square meter (wattsm).
        // average daily global radiation, W m-2.
        let wattsm = GetFromClim(CLIMATE_METRIC_IRRD, Daynum) * 697.45 / (DayLength * 60.);
        // Compute pstand as an empirical function of wattsm (based on Baker et al., 1972).
        // gross photosynthesis for a non-stressed full canopy.
        let pstand = 2.3908 + wattsm * (1.37379 - wattsm * 0.00054136);
        // Convert it to gross photosynthesis per plant (pplant), using PerPlantArea and corrections for light interception by canopy, ambient CO2 concentration, water stress and low N in the leaves.
        let mut pplant = 0.;
        // actual gross photosynthetic rate, g per plant per day.
        match profile.light_intercept_method {
            LightInterceptMethod::Latered => {
                let mut pstand_remain = pstand;
                for i in (0..20).rev() {
                    if pstand_remain <= 0. {
                        break;
                    }
                    if LightInterceptLayer[i] <= 0. {
                        continue;
                    }
                    let page = 1. - (AverageLeafAge[i] / drop_leaf_age(LeafArea[i])).powi(2);
                    let mut pplant_inc = 0.001
                        * pstand_remain
                        * LightInterceptLayer[i]
                        * PerPlantArea
                        * profile.ptsred
                        * pnetcor
                        * ptnfac
                        * page;
                    if pplant_inc > pstand_remain {
                        pplant_inc = pstand_remain;
                    }
                    pplant += pplant_inc;
                    pstand_remain -= pplant_inc;
                }
            }
            _ => {
                pplant = 0.001
                    * pstand
                    * LightIntercept
                    * PerPlantArea
                    * profile.ptsred
                    * pnetcor
                    * ptnfac;
            }
        }
        // Compute the photorespiration factor (rsubl) as a linear function af average day time temperature.
        let rsubl = 0.0032125 + 0.0066875 * DayTimeTemp; // photorespiration factor.

        // Photorespiration (lytres) is computed as a proportion of gross photosynthetic rate.
        let lytres = rsubl * pplant; // rate of photorespiration, g per plant per day.

        // Old stems are those more than voldstm = 32 calendar days old.
        // Maintenance respiration is computed on the basis of plant dry weight, minus the old stems and the dry tissue of opened bolls.
        let voldstm = 32;
        let kkday = Kday - voldstm; // day of least recent actively growing stems.

        // weight of old stems.
        let oldstmwt = if kkday < 1 {
            0.
        } else {
            StemWeight[kkday as usize]
        };
        // maintenance respiration, g per plant per day.
        let bmain = (PlantWeight - CottonWeightOpenBolls - BurrWeightOpenBolls - oldstmwt) * rsubo;
        // Net photosynthesis is computed by substracting photo-respiration and maintenance respiration from the gross rate of photosynthesis. To avoid computational problems, make sure that pts is positive and non-zero.
        // intermediate computation of NetPhotosynthesis.
        let mut pts = pplant - lytres - bmain;
        if pts < 0.00001 {
            pts = 0.00001;
        }
        // The growth respiration (gsubr) supplies energy for converting the supplied carbohydrates to plant tissue dry
        // matter. 0.68182 converts CO2 to CH2O. NetPhotosynthesis is the computed net photosynthesis, in g per plant
        // per day.

        NetPhotosynthesis = pts / (1. + gsubr) * 0.68182;
        // CumNetPhotosynth is the cumulative value of NetPhotosynthesis, from day of emergence.
        CumNetPhotosynth += NetPhotosynthesis;
        Ok(())
    }
    /// This function sets the water saturation of the soil layers below the water table, if it has been defined in the input.
    /// It is called from [Profile::soil_procedures()] if water table data have been input.
    ///
    /// The following global variables are referenced here:
    /// * [Daynum]
    /// * [dl]
    /// * [MaxWaterCapacity]
    /// * [nk]
    /// * [nl]
    /// * [PoreSpace]
    /// * [MaxWaterCapacity]
    /// * [RowSpace]
    /// * [wk]
    ///
    /// The following global variables are set here:
    /// * [addwtbl]
    /// * [ElCondSatSoilToday]
    /// * [WaterTableLayer]
    /// * [VolWaterContent]
    fn watertable(&mut self, profile: &Profile) -> Result<(), Cotton2KError> {
        if profile.num_watertable_data == 0 {
            return Ok(());
        }
        // Find the depth of water table for this day.
        let mut lwtable = 201f64; // level of water table on this day, cm
        for ao in &profile.agronomy_operations {
            match ao {
                AgronomyOperation::watertable { date, level, ecs } => unsafe {
                    if Daynum as u32 >= date.ordinal() {
                        lwtable = *level;
                        ElCondSatSoilToday = *ecs;
                    }
                },
                _ => {}
            }
        }
        // Find the number of the uppermost layer of water table
        unsafe {
            WaterTableLayer = 1000;
            let mut sumdl = 0f64; // sum of depth of consecutive soil layers
            for l in 0..nl as usize {
                sumdl += dl[l];
                if sumdl > lwtable {
                    WaterTableLayer = l as i32;
                    break;
                }
            }
        }
        // The total water entering the soil slab (addwtbl) is computed. It is used to check the water balance in the soil.

        unsafe {
            for l in 0..nl as usize {
                if l as i32 >= WaterTableLayer {
                    for k in 0..nk as usize {
                        // previous water content of a cell
                        let vh2ocx = VolWaterContent[l][k];
                        VolWaterContent[l][k] = PoreSpace[l];
                        addwtbl +=
                            10. * (VolWaterContent[l][k] - vh2ocx) * dl[l] * wk[k] / RowSpace;
                    }
                } else {
                    // Make sure that (in case water table was lowered) water content is not higher than MaxWaterCapacity and adjust addwtbl.
                    for k in 0..nk as usize {
                        if VolWaterContent[l][k] > MaxWaterCapacity[l] {
                            // previous water content of a cell
                            let vh2ocx = VolWaterContent[l][k];
                            VolWaterContent[l][k] = MaxWaterCapacity[l];
                            addwtbl +=
                                10. * (VolWaterContent[l][k] - vh2ocx) * dl[l] * wk[k] / RowSpace;
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

/// This function computes the water redistribution in the soil or surface irrigation (by flooding or sprinklers).
/// It is called by [SoilProcedures()].
/// It calls function [Drain()].
///
/// The following argument is used:
/// - applyWat = amount of water applied, mm.
///
/// The following global variables are referenced:
/// * [dl]
/// * [nk]
/// * [RowSpace].
///
/// The following global variables are set:
/// * [CumWaterDrained]
/// * [VolWaterContent]
unsafe fn GravityFlow(applywat: f64) {
    // Add the applied amount of water to the top soil cell of each column.
    for k in 0..nk as usize {
        VolWaterContent[0][k] += 0.10 * applywat / dl[0];
    }
    // Call function Drain() to compute downflow of water.
    // water drained out of the slab, mm.
    let WaterDrainedOut: f64 = Drain();
    // If there is drainage out of the slab, transform it to mm, and update the cumulative drainage (CumWaterDrained)
    if WaterDrainedOut > 0. {
        CumWaterDrained += 10. * WaterDrainedOut / RowSpace;
    }
}

/// This function simulates the leaf water potential of cotton plants. It has been adapted from the model of Moshe Meron (The relation of cotton leaf water potential to soil water content in the irrigated management range. PhD dissertation, UC Davis, 1984).
///
///  It is called from [Profile::stress()]. It calls [wcond()] and [LeafResistance()].
///
/// The following global variables are referenced here:
/// * [AgeOfPreFruNode]
/// * [AverageSoilPsi]
/// * [beta]
/// * [dl]
/// * [Kday]
/// * [LeafAge]
/// * [NumFruitBranches]
/// * [NumLayersWithRoots]
/// * [NumNodes]
/// * [NumPreFruNodes]
/// * [NumVegBranches]
/// * [PlantHeight]
/// * [PoreSpace]
/// * [ReferenceETP]
/// * [RootColNumLeft]
/// * [RootColNumRight]
/// * [RootWtCapblUptake]
/// * [SaturatedHydCond]
/// * [SoilPsi]
/// * [thad]
/// * [thts]
/// * [VolWaterContent]
/// * [wk]
///
/// The following global variables are set here:
/// * [LwpMin]
/// * [LwpMax]
fn LeafWaterPotential() {
    // Constant parameters used:
    const cmg: f64 = 3200.; // length in cm per g dry weight of roots, based on an average root diameter of 0.06 cm, and a specific weight of 0.11 g  dw per cubic cm.
    const psild0: f64 = -1.32; // maximum values of LwpMin
    const psiln0: f64 = -0.40; // maximum values of LwpMax.
    const rtdiam: f64 = 0.06; // average root diameter in cm.
    const vpsil: [f64; 13] = [
        0.48, -5.0, 27000., 4000., 9200., 920., 0.000012, -0.15, -1.70, -3.5, 0.1e-9, 0.025, 0.80,
    ];
    // Leaf water potential is not computed during 10 days after emergence. Constant values are assumed for this period.
    unsafe {
        if Kday <= 10 {
            LwpMax = psiln0;
            LwpMin = psild0;
            return;
        }
    }
    // Compute shoot resistance (rshoot) as a function of plant height.
    let rshoot: f64 = vpsil[0] * unsafe { PlantHeight } / 100.; // shoot resistance, Mpa hours per cm.

    // Assign zero to summation variables
    let mut psinum = 0f64; // sum of RootWtCapblUptake for all soil cells with roots.
    let mut rootvol = 0f64; // sum of volume of all soil cells with roots.
    let mut rrlsum = 0f64; // weighted sum of reciprocals of rrl.
    let rroot; // root resistance, Mpa hours per cm.
    let mut sumlv = 0f64; // weighted sum of root length, cm, for all soil cells with roots.
    let mut vh2sum = 0f64; // weighted sum of soil water content, for all soil cells with roots.

    // Loop over all soil cells with roots. Check if RootWtCapblUptake is greater than vpsil[10].
    // All average values computed for the root zone, are weighted by RootWtCapblUptake (root weight capable of uptake), but the weight assigned will not be greater than vpsil[11].
    unsafe {
        for l in 0..NumLayersWithRoots as usize {
            for k in (RootColNumLeft[l] as usize)..RootColNumLeft[l] as usize {
                if RootWtCapblUptake[l][k] >= vpsil[10] {
                    psinum += fmin(RootWtCapblUptake[l][k], vpsil[11]);
                    sumlv += fmin(RootWtCapblUptake[l][k], vpsil[11]) * cmg;
                    rootvol += dl[l] * wk[k];
                    // root resistance per g of active roots.
                    let rrl = if SoilPsi[l][k] <= vpsil[1] {
                        vpsil[2] / cmg
                    } else {
                        (vpsil[3] - SoilPsi[l][k] * (vpsil[4] + vpsil[5] * SoilPsi[l][k])) / cmg
                    };
                    rrlsum += fmin(RootWtCapblUptake[l][k], vpsil[11]) / rrl;
                    vh2sum += VolWaterContent[l][k] * fmin(RootWtCapblUptake[l][k], vpsil[11]);
                }
            }
        }
    }
    // Compute average root resistance (rroot) and average soil water content (vh2).
    let dumyrs: f64; // intermediate variable for computing cond.
    let vh2: f64; // average of soil water content, for all soil soil cells with roots.
    if psinum > 0. && sumlv > 0. {
        rroot = psinum / rrlsum;
        vh2 = vh2sum / psinum;
        dumyrs = fmax(
            1.001,
            (1. / (std::f64::consts::PI * sumlv / rootvol)).sqrt() / rtdiam,
        );
    } else {
        rroot = 0.;
        vh2 = unsafe { thad[0] };
        dumyrs = 1.001;
    }
    // Compute hydraulic conductivity (cond), and soil resistance near the root surface (rsoil). soil hydraulic conductivity near the root surface.
    let cond = fmax(
        unsafe {
            wcond(
                vh2,
                thad[0],
                thts[0],
                beta[0],
                SaturatedHydCond[0],
                PoreSpace[0],
            )
        } / 24.
            * 2.
            * sumlv
            / rootvol
            / dumyrs.ln(),
        vpsil[6],
    );
    let rsoil = 0.0001 / (2. * std::f64::consts::PI * cond); // soil resistance, Mpa hours per cm.

    // Compute leaf resistance (LeafResistance) as the average of the resistances of all existing leaves. The resistance of an individual leaf is a function of its age. Function LeafResistance is called to compute it.
    // This is executed for all the leaves of the plant.
    let mut numl = 0; // number of leaves.
    let mut sumrl = 0f64; // sum of leaf resistances for all the plant.
    unsafe {
        for j in 0..NumPreFruNodes as usize
        // loop prefruiting nodes
        {
            numl += 1;
            sumrl += LeafResistance(AgeOfPreFruNode[j]);
        }
    }
    //
    let mut nbrch: i32; // number of fruiting branches on a vegetative branch.
    let mut nnid: i32; // number of nodes on a fruiting branch.
    unsafe {
        for k in 0..NumVegBranches as usize {
            // loop for all other nodes
            nbrch = NumFruitBranches[k];
            for l in 0..nbrch as usize {
                nnid = NumNodes[k][l];
                for m in 0..nnid as usize {
                    numl += 1;
                    sumrl += LeafResistance(LeafAge[k][l][m]);
                }
            }
        }
    }
    let rleaf = sumrl / numl as f64; // leaf resistance, Mpa hours per cm.

    // The total resistance to transpiration, MPa hours per cm, (rtotal) is computed.
    let rtotal = rsoil + rroot + rshoot + rleaf;
    // Compute maximum (early morning) leaf water potential, LwpMax, from soil water potential (AverageSoilPsi, converted from bars to MPa).
    // Check for minimum and maximum values.
    unsafe {
        LwpMax = vpsil[7] + 0.1 * AverageSoilPsi;
        if LwpMax < vpsil[8] {
            LwpMax = vpsil[8];
        }
        if LwpMax > psiln0 {
            LwpMax = psiln0;
        }
    }
    //     Compute minimum (at time of maximum transpiration rate) leaf water
    //     potential, LwpMin, from
    //  maximum transpiration rate (etmax) and total resistance to transpiration
    //  (rtotal).
    let mut etmax = 0f64; // the maximum hourly rate of evapotranspiration for this day.
    for ihr in 0..24 {
        //  hourly loop
        unsafe {
            if ReferenceETP[ihr] > etmax {
                etmax = ReferenceETP[ihr];
            }
        }
    }
    unsafe {
        LwpMin = LwpMax - 0.1 * fmax(etmax, vpsil[12]) * rtotal;
        //     Check for minimum and maximum values.
        if LwpMin < vpsil[9] {
            LwpMin = vpsil[9];
        }
        if LwpMin > psild0 {
            LwpMin = psild0;
        }
    }
}
/// This function computes the water redistribution in the soil after irrigation by a drip system.
/// It also computes the resulting redistribution of nitrate and urea N.
/// It is called by SoilProcedures() noitr times per day.
/// It calls function CellDistrib().
///
/// The following argument is used:
/// Drip - amount of irrigation applied by the drip method, mm.
///
/// The following global variables are referenced:
/// * [dl]
/// * [LocationColumnDrip]
/// * [LocationLayerDrip]
/// * [MaxWaterCapacity]
/// * [nk]
/// * [nl]
/// * [NO3FlowFraction]
/// * [PoreSpace]
/// * [RowSpace]
/// * [wk]
///
/// The following global variables are set:
/// * [CumWaterDrained]
/// * [SoilNitrogenLoss]
/// * [VolWaterContent]
/// * [VolNo3NContent]
/// * [VolUreaNContent]
fn DripFlow(Drip: f64) -> Result<(), Cotton2KError> {
    let mut dripw: [f64; 40] = [0f64; 40]; // amount of water applied, or going from one ring of
                                           // soil cells to the next one, cm3. (array)
    let mut dripn: [f64; 40] = [0f64; 40]; // amount of nitrate N applied, or going from one ring of
                                           // soil soil cells to the next one, mg. (array)
    let mut dripu: [f64; 40] = [0f64; 40]; // amount of urea N applied, or going from one ring of
                                           // soil soil cells to the next one, mg. (array)

    // Incoming flow of water (Drip, in mm) is converted to dripw(0), in cm^3 per slab.
    dripw[0] = Drip * unsafe { RowSpace } * 0.10;
    // Wetting the cell in which the emitter is located.
    let mut h2odef; // the difference between the maximum water capacity (at a water content of uplimit) of this ring of soil cell, and the actual water content, cm3.
    let l0 = unsafe { LocationLayerDrip } as usize; //  layer where the drip emitter is situated
    let k0 = unsafe { LocationColumnDrip } as usize; //  column where the drip emitter is situated

    // It is assumed that wetting cannot exceed MaxWaterCapacity of this cell.
    // Compute h2odef, the amount of water needed to saturate this cell.
    h2odef = unsafe { (MaxWaterCapacity[l0] - VolWaterContent[l0][k0]) * dl[l0] * wk[k0] };
    // If maximum water capacity is not exceeded - update VolWaterContent of this cell and exit the function.
    if dripw[0] <= h2odef {
        unsafe {
            VolWaterContent[l0][k0] += dripw[0] / (dl[l0] * wk[k0]);
        }
        return Ok(());
    }
    // If maximum water capacity is exceeded - calculate the excess of water flowing out of this cell (in cm3 per slab) as dripw[1].
    // The next ring of cells (kr=1) will receive it as incoming water flow.
    dripw[1] = dripw[0] - h2odef;
    // Compute the movement of nitrate N to the next ring
    let mut cnw; // concentration of nitrate N in the outflowing water
    if unsafe { VolNo3NContent[l0][k0] } > 1.0e-30 {
        cnw = unsafe {
            VolNo3NContent[l0][k0] / (VolWaterContent[l0][k0] + dripw[0] / (dl[l0] * wk[k0]))
        };
        // cnw is multiplied by dripw[1] to get dripn[1], the amount of nitrate N going out to the next ring of cells.
        // It is assumed, however, that not more than a proportion (NO3FlowFraction) of the nitrate N in this cell can be removed in one iteration.
        if unsafe { (cnw * MaxWaterCapacity[l0]) < (NO3FlowFraction[l0] * VolNo3NContent[l0][k0]) }
        {
            dripn[1] = unsafe { NO3FlowFraction[l0] * VolNo3NContent[l0][k0] * dl[l0] * wk[k0] };
            unsafe { VolNo3NContent[l0][k0] = (1. - NO3FlowFraction[l0]) * VolNo3NContent[l0][k0] };
        } else {
            dripn[1] = dripw[1] * cnw;
            unsafe {
                VolNo3NContent[l0][k0] = MaxWaterCapacity[l0] * cnw;
            }
        }
    }
    // The movement of urea N to the next ring is computed similarly.
    let mut cuw; // concentration of urea N in the outflowing water
    if unsafe { VolUreaNContent[l0][k0] } > 1.0e-30 {
        cuw = unsafe {
            VolUreaNContent[l0][k0] / (VolWaterContent[l0][k0] + dripw[0] / (dl[l0] * wk[k0]))
        };
        if unsafe { (cuw * MaxWaterCapacity[l0]) < (NO3FlowFraction[l0] * VolUreaNContent[l0][k0]) }
        {
            dripu[1] = unsafe { NO3FlowFraction[l0] * VolUreaNContent[l0][k0] * dl[l0] * wk[k0] };
            unsafe {
                VolUreaNContent[l0][k0] = (1. - NO3FlowFraction[l0]) * VolUreaNContent[l0][k0];
            }
        } else {
            dripu[1] = dripw[1] * cuw;
            unsafe {
                VolUreaNContent[l0][k0] = MaxWaterCapacity[l0] * cuw;
            }
        }
    }
    let mut defcit = ndarray::Array2::<f64>::zeros((40, 20)); // array of the difference between water capacity and actual water content in each cell of the ring

    // Set VolWaterContent of the cell in which the drip is located to MaxWaterCapacity.
    unsafe {
        VolWaterContent[l0][k0] = MaxWaterCapacity[l0];
    }
    // Loop of concentric rings of cells, starting from ring 1.
    // Assign zero to the sums sv, st, sn, sn1, su and su1.
    for kr in 1..maxl as usize {
        let mut uplimit; //  upper limit of soil water content in a soil cell
        let mut sv = 0f64; // sum of actual water content in a ring of cells, cm3
        let mut st = 0f64; // sum of total water capacity in a ring of cells, cm3
        let mut sn = 0f64; // sum of nitrate N content in a ring of cells, mg.
        let mut sn1 = 0f64; // sum of movable nitrate N content in a ring of cells, mg
        let mut su = 0f64; // sum of urea N content in a ring of cells, mg
        let mut su1 = 0f64; // sum of movable urea N content in a ring of cells, mg
        let radius = (6 * kr) as f64; // radius (cm) of the wetting ring
        let mut dist; // distance (cm) of a cell center from drip location

        for l in 1..unsafe { nl as usize } {
            // Upper limit of water content is the porespace volume in layers below the water table, MaxWaterCapacity in other layers.
            if l >= unsafe { WaterTableLayer as usize } {
                uplimit = unsafe { PoreSpace[l] };
            } else {
                unsafe {
                    uplimit = MaxWaterCapacity[l];
                }
            }
            for k in 0..unsafe { nk as usize } {
                // Compute the sums sv, st, sn, sn1, su and su1 within the radius limits of this ring.
                // The array defcit is the sum of difference between uplimit and VolWaterContent of each cell.
                dist = cell_distance(l, k, l0, k0, unsafe { RowSpace })?;
                if dist <= radius && dist > (radius - 6.) {
                    sv += unsafe { VolWaterContent[l][k] * dl[l] * wk[k] };
                    st += uplimit * unsafe { dl[l] * wk[k] };
                    sn += unsafe { VolNo3NContent[l][k] * dl[l] * wk[k] };
                    sn1 += unsafe { VolNo3NContent[l][k] * dl[l] * wk[k] * NO3FlowFraction[l] };
                    su += unsafe { VolUreaNContent[l][k] * dl[l] * wk[k] };
                    su1 += unsafe { VolUreaNContent[l][k] * dl[l] * wk[k] * NO3FlowFraction[l] };
                    defcit[[l, k]] = uplimit - unsafe { VolWaterContent[l][k] };
                } else {
                    defcit[[l, k]] = 0.;
                }
            }
        }
        // Compute the amount of water needed to saturate all the cells in this ring (h2odef).
        h2odef = st - sv;
        // Test if the amount of incoming flow, dripw(kr), is greater than h2odef.
        if dripw[kr] <= h2odef {
            // In this case, this will be the last wetted ring.
            // Update VolWaterContent in this ring, by wetting each cell in proportion to its defcit.
            // Update VolNo3NContent and VolUreaNContent of the cells in this ring by the same proportion.
            // This is executed for all the cells in the ring.
            for l in 1..unsafe { nl } as usize {
                for k in 0..unsafe { nk } as usize {
                    dist = cell_distance(l, k, l0, k0, unsafe { RowSpace })?;
                    if dist <= radius && dist > (radius - 6.) {
                        unsafe {
                            VolWaterContent[l][k] += dripw[kr] * defcit[[l, k]] / h2odef;
                            VolNo3NContent[l][k] += dripn[kr] * defcit[[l, k]] / h2odef;
                            VolUreaNContent[l][k] += dripu[kr] * defcit[[l, k]] / h2odef;
                        }
                    }
                }
            }
            return Ok(());
        }
        // If dripw(kr) is greater than h2odef, calculate cnw and cuw as the concentration of nitrate and urea N in the total water of this ring after receiving the incoming water and nitrogen.
        cnw = (sn + dripn[kr]) / (sv + dripw[kr]);
        cuw = (su + dripu[kr]) / (sv + dripw[kr]);
        let drwout = dripw[kr] - h2odef; //  the amount of water going out of a ring, cm3.

        // Compute the nitrate and urea N going out of this ring, and their amount lost from this ring.
        // It is assumed that not more than a certain part of the total nitrate or urea N (previously computed as sn1 an su1) can be lost from a ring in one iteration.
        // drnout and xnloss are adjusted accordingly.
        // druout and xuloss are computed similarly for urea N.
        let mut drnout = drwout * cnw; //  the amount of nitrate N going out of a ring, mg
        let mut xnloss = 0f64; // the amount of nitrate N lost from a ring, mg
        if drnout > dripn[kr] {
            xnloss = drnout - dripn[kr];
            if xnloss > sn1 {
                xnloss = sn1;
                drnout = dripn[kr] + xnloss;
            }
        }
        let mut druout = drwout * cuw; //  the amount of urea N going out of a ring, mg
        let mut xuloss = 0f64; // the amount of urea N lost from a ring, mg
        if druout > dripu[kr] {
            xuloss = druout - dripu[kr];
            if xuloss > su1 {
                xuloss = su1;
                druout = dripu[kr] + xuloss;
            }
        }
        // For all the cells in the ring, as in the 1st cell, saturate VolWaterContent to uplimit, and update VolNo3NContent and VolUreaNContent.
        for l in 1..unsafe { nl } as usize {
            uplimit = unsafe {
                if l >= WaterTableLayer as usize {
                    PoreSpace[l]
                } else {
                    MaxWaterCapacity[l]
                }
            };
            for k in 0..unsafe { nk } as usize {
                dist = cell_distance(l, k, l0, k0, unsafe { RowSpace })?;
                if dist <= radius && dist > (radius - 6.) {
                    unsafe {
                        VolWaterContent[l][k] = uplimit;
                        VolNo3NContent[l][k] = if xnloss <= 0. {
                            uplimit * cnw
                        } else {
                            VolNo3NContent[l][k] * (1. - xnloss / sn)
                        };
                        VolUreaNContent[l][k] = if xuloss <= 0. {
                            uplimit * cuw
                        } else {
                            VolUreaNContent[l][k] * (1. - xuloss / su)
                        };
                    }
                }
            }
        }
        // The outflow of water, nitrate and urea from this ring will be the inflow into the next ring.
        if kr < (unsafe { nl } as usize - l0 - 1) && kr < maxl as usize - 1 {
            dripw[kr + 1] = drwout;
            dripn[kr + 1] = drnout;
            dripu[kr + 1] = druout;
        } else
        // If this is the last ring, the outflowing water will be added to the drainage, CumWaterDrained, the outflowing nitrogen to SoilNitrogenLoss.
        {
            unsafe {
                CumWaterDrained += 10. * drwout / RowSpace;
                SoilNitrogenLoss += drnout + druout;
            }
            return Ok(());
        }
        // Repeat all these procedures for the next ring.
    }
    Ok(())
}

fn drop_leaf_age(lai: f64) -> f64 {
    140. - 1. * lai
}
