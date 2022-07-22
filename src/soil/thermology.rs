use crate::atmosphere::{clearskyemiss, vapor_pressure};
use crate::{
    albedo, bEnd, dl, es1hour, es2hour, isw, nk, nl, thad, wk, ActualSoilEvaporation,
    ActualTranspiration, AirTemp, CanopyBalance, Clim, CloudCoverRatio, CloudTypeCorr,
    Cotton2KError, CumEvaporation, DayEndMulch, DayOfSimulation, DayPlant, DayStart, DayStartMulch,
    Daynum, DeepSoilTemperature, FieldCapacity, FoliageTemp, HeatCapacitySoilSolid, Kday,
    MulchIndicator, MulchSurfaceBalance, MulchTemp, MulchTranLW, MulchTranSW, PlantHeight,
    PlantRowColumn, PoreSpace, PredictEmergence, Radiation, ReferenceETP, ReferenceTransp,
    RelativeHumidity, RowSpace, Scratch21, SensibleHeatTransfer, SitePar, SoilSurfaceBalance,
    SoilTemp, SoilTempDailyAvrg, ThermalCondSoil, VolWaterContent, WindSpeed,
};
const maxl: usize = 40;

#[derive(Debug, Clone, Copy)]
pub struct SoilThermology {
    /// the relative radiation received by a soil column, as affected by shading by plant canopy.
    pub rracol: [f64; 20],
    numiter: u64,
    /// equal to the dl array in a columnn, or wk in a row.
    dz: [f64; maxl],
    /// array of previous soil temperatures.
    ts0: [f64; maxl],
    /// array of soil temperatures.
    ts1: [f64; maxl],
    /// heat capacity of soil layer (cal cm-3 oC-1).
    hcap: [f64; maxl],
}

impl SoilThermology {
    /// Called from soil_thermology() at the start of the simulation.
    /// It sets initial values to soil and canopy temperatures.
    ///
    /// The following global variables are referenced here:
    /// Clim (structure), DayFinish, Daynum, DayStart, nl, SitePar.
    ///
    /// The following global variables are set here:
    /// DeepSoilTemperature, SoilTemp.
    pub fn new() -> Self {
        // If there is an output flag for soil temperatures, an error message pops up for defining the start and stop dates for this output.
        // Compute initial values of soil temperature: It is assumed that at the start of simulation the temperature of the first soil layer (upper boundary) is equal to the average air temperature of the previous five days (if climate data not available - start from first climate data).
        //
        // NOTE: For a good simulation of soil temperature, it is recommended to start simulation at least 10 days before planting date.
        // This means that climate data should be available for this period. This is especially important if emergence date has to be simulated.
        let mut idd = unsafe { Daynum - 4 - DayStart }; // number of days minus 4 from start of simulation.
        if idd < 0 {
            idd = 0;
        }
        let mut tsi1 = 0.; // Upper boundary (surface layer) initial soil temperature, C.
        for i in idd as usize..(idd + 5) as usize {
            tsi1 += unsafe { Clim[i].Tmax + Clim[i].Tmin };
        }
        tsi1 /= 10.;
        // The temperature of the last soil layer (lower boundary) is computed as a sinusoidal function of day of year, with site-specific parameters.
        unsafe {
            DeepSoilTemperature = SitePar[9]
                + SitePar[10]
                    * (2. * std::f64::consts::PI * (Daynum as f64 - SitePar[11]) / 365.).sin();
        }
        // SoilTemp is assigned to all columns, converted to degrees K.
        tsi1 += 273.161;
        unsafe {
            DeepSoilTemperature += 273.161;
        }
        for l in 0..unsafe { nl } as usize {
            // The temperatures of the other soil layers are linearly interpolated.
            // computed initial soil temperature, C, for each layer
            let tsi = ((unsafe { nl } as usize - l - 1) as f64 * tsi1
                + l as f64 * unsafe { DeepSoilTemperature })
                / (unsafe { nl } - 1) as f64;
            for k in 0..unsafe { nk } as usize {
                unsafe {
                    SoilTemp[l][k] = tsi;
                }
            }
        }
        SoilThermology {
            rracol: [1.; 20],
            numiter: 0,
            dz: [0.; 40],
            ts0: [0.; 40],
            ts1: [0.; 40],
            hcap: [0.; 40],
        }
    }

    /// This is the main part of the soil temperature sub-model. It is called daily from SimulateThisDay(). It calls the following functions:
    /// EnergyBalance(), PredictEmergence(), SoilHeatFlux(), SoilTemperatureInit().
    ///
    /// The following global variables are referenced here:
    /// ActualTranspiration, AirTemp, Date, DayEndMulchDaynum, DayPlant,
    /// DayStart, DayStartMulch, DewPointTemp, dl, FoliageTemp, isw,
    /// MulchIndicator, MulchTemp, nk, nl, PlantRowColumn,
    /// ReferenceETP, ReferenceTransp, RelativeHumidity, RowSpace,
    /// rracol, SitePar, thad, wk
    ///
    /// The following global variables are set here:
    /// ActualSoilEvaporation, bEnd, CumEvaporation, DeepSoilTemperature, es,
    /// SoilTemp, SoilTempDailyAvrg, VolWaterContent,
    pub unsafe fn soil_thermology(&mut self) {
        //     Compute dts, the daily change in deep soil temperature (C), as
        //  a site-dependent function of Daynum.
        let dts = 2. * std::f64::consts::PI * SitePar[10] / 365.
            * (2. * std::f64::consts::PI * (Daynum as f64 - SitePar[11]) / 365.).cos();
        //     Define iter1 and dlt for hourly time step.
        let iter1 = 24; // number of iterations per day.
        let dlt = 3600.; // time (seconds) of one iteration.
        let mut kk = 1; // number of soil columns for executing computations.
                        //     If there is no canopy cover, no horizontal heat flux is assumed, kk
                        //     = 1.
                        //  Otherwise it is equal to the number of columns in the slab.
        let mut shadeav = 0.; // average shaded area in all shaded soil columns.
                              //     isw defines the type of soil temperature computation.
        if isw > 1 {
            let mut shadetot = 0.; // sum of shaded area in all shaded soil columns.
            let mut nshadedcol = 0; // number of at least partially shaded soil columns.
            kk = nk as usize;
            for k in 0..nk as usize {
                if self.rracol[k] <= 0.99 {
                    shadetot += 1. - self.rracol[k];
                    nshadedcol += 1;
                }
            }
            if nshadedcol > 0 {
                shadeav = shadetot / nshadedcol as f64;
            }
        }
        //     Set daily averages of soil temperature to zero.
        for l in 0..nl as usize {
            for k in 0..nk as usize {
                SoilTempDailyAvrg[l][k] = 0.;
            }
        }
        //     es and ActualSoilEvaporation are computed as the average for the
        //     whole soil
        //  slab, weighted by column widths.
        let mut es = 0f64; // potential evaporation rate, mm day-1
        ActualSoilEvaporation = 0.;
        //     Start hourly loop of iterations.
        for ihr in 0..iter1 {
            //     Update the temperature of the last soil layer (lower boundary
            //     conditions).
            DeepSoilTemperature += dts * dlt / 86400.;
            // actual transpiration (mm s-1) for this hour
            let etp0 = if ReferenceTransp > 0.000001 {
                ActualTranspiration * ReferenceETP[ihr] / ReferenceTransp / dlt
            } else {
                0.
            };

            // Compute vertical transport for each column
            for k in 0..kk as usize {
                //     Set SoilTemp for the lowest soil layer.
                SoilTemp[(nl - 1) as usize][k] = DeepSoilTemperature;
                //     Compute transpiration from each column, weighted by its
                //     relative shading.

                // actual hourly transpiration (mm s-1) for a column.
                let etp1 = if shadeav > 0.000001 {
                    etp0 * (1. - self.rracol[k]) / shadeav
                } else {
                    0.
                };
                //     Check if mulch is on for this date and for this column.
                // is true if this column is covered with plastic mulch now, false if not.
                let bMulchon =
                    if MulchIndicator == 0 || Daynum < DayStartMulch || Daynum > DayEndMulch {
                        false
                    } else {
                        if MulchIndicator == 1 {
                            true
                        } else if MulchIndicator == 2 {
                            if k >= PlantRowColumn as usize && k <= (PlantRowColumn + 1) as usize {
                                false
                            } else {
                                true
                            }
                        } else if MulchIndicator >= 3 {
                            if k >= (PlantRowColumn - 1) as usize
                                && k <= (PlantRowColumn + 2) as usize
                            {
                                false
                            } else {
                                true
                            }
                        } else {
                            false
                        }
                    };
                let mut ess = 0.; //   evaporation rate from surface of a soil column (mm / sec).
                if !bMulchon {
                    // The potential evaporation rate (escol1k) from a column is the sum of the radiation
                    // component of the Penman equation(es1hour), multiplied by the
                    // relative radiation reaching this column, and the wind and
                    // vapor deficit component of the Penman equation (es2hour).
                    // potential evaporation fron soil surface of a column, mm per hour.
                    let mut escol1k = es1hour[ihr] * self.rracol[k] + es2hour[ihr];
                    es += escol1k * wk[k];
                    // Compute actual evaporation from soil surface. Update VolWaterContent of the soil soil cell, and add to daily sum of actual evaporation.
                    let evapmax = 0.9 * (VolWaterContent[0][k] - thad[0]) * 10. * dl[0]; // maximum possible evaporatio from a soil cell near the surface.
                    if escol1k > evapmax {
                        escol1k = evapmax;
                    }
                    VolWaterContent[0][k] -= 0.1 * escol1k / dl[0];
                    ActualSoilEvaporation += escol1k * wk[k];
                    ess = escol1k / dlt;
                }
                // Call EnergyBalance to compute soil surface and canopy temperature.
                self.soil_energy_balance(ihr, k, bMulchon, ess, etp1);
                if bEnd {
                    return;
                }
            }

            // Compute soil temperature flux in the vertical direction.
            // Assign iv = 1, layer = 0, nn = nl.
            let mut iv = 1; // indicates vertical (=1) or horizontal (=0) flux.
            let mut nn = nl as usize; // number of array members for heat flux.
            let mut layer = 0; // soil layer number
            let mut tsolav: [f64; 40] = [0.; 40]; // hourly average soil temperature C, of a soil layer.
            for k in 0..kk {
                // Loop over kk columns, and call SoilHeatFlux().
                self.heat_flux(dlt, iv, nn, layer, k);
            }
            // If no horizontal heat flux is assumed, make all array members of SoilTemp equal to the value computed for the first column. Also, do the same for array memebers of VolWaterContent.
            if isw <= 1 {
                for l in 0..nl as usize {
                    for k in 1..nk as usize {
                        SoilTemp[l][k] = SoilTemp[l][0];
                        if l == 0 {
                            VolWaterContent[l][k] = VolWaterContent[l][0];
                        }
                    }
                }
            }

            // Compute horizontal transport for each layer
            //
            // Compute soil temperature flux in the horizontal direction, when
            // isw = 2. Assign iv = 0 and nn = nk. Start loop for soil layers,
            // and call SoilHeatFlux.
            if isw > 1 {
                iv = 0;
                nn = nk as usize;
                for l in 0..nl as usize {
                    layer = l;
                    self.heat_flux(dlt, iv, nn, layer, l);
                }
            }
            //     Compute average temperature of soil layers, in degrees C.
            for l in 0..nl as usize {
                for k in 0..nk as usize {
                    SoilTempDailyAvrg[l][k] += SoilTemp[l][k];
                    tsolav[l] += SoilTemp[l][k] - 273.161;
                }
                tsolav[l] = tsolav[l] / nk as f64;
            }
            // If emergence date is to be simulated, call PredictEmergence().
            if isw == 0 && Daynum >= DayPlant {
                PredictEmergence(ihr as i32);
            }
        }
        // At the end of the day compute actual daily evaporation and its cumulative sum.
        if kk == 1 {
            es /= wk[1];
            ActualSoilEvaporation = ActualSoilEvaporation / wk[1];
        } else {
            es /= RowSpace;
            ActualSoilEvaporation = ActualSoilEvaporation / RowSpace;
        }
        CumEvaporation += ActualSoilEvaporation;
        if Kday > 0 {
            Scratch21[(DayOfSimulation - 1) as usize].es = es;
            Scratch21[(DayOfSimulation - 1) as usize].cumEvaporation = CumEvaporation;
        }
        // compute daily averages.
        for l in 0..nl as usize {
            for k in 0..nk as usize {
                SoilTempDailyAvrg[l][k] = SoilTempDailyAvrg[l][k] / iter1 as f64;
            }
        }
    }
    /// Solves the energy balance equations at the soil surface, and at the foliage/atmosphere interface.
    /// It computes the resulting temperatures of the soil surface, plastic mulch (if exists) and the plant canopy.
    ///
    /// Units for all energy fluxes are: cal cm-2 sec-1.
    ///
    /// It is called from soil_thermology(), on each hourly time step and for each soil column.
    /// It calls functions [clearskyemiss()], [vapor_pressure()], SensibleHeatTransfer(), SoilMulchBalance(), SoilSurfaceBalance() and CanopyBalance().
    ///
    /// The following arguments are used in this function:
    /// * `ihr` - the time of day in hours.
    /// * `k` - soil column number.
    /// * `bMulchon` - is true if this column is covered with plastic mulch, false if not.
    /// * `ess` - evaporation from surface of a soil column (mm / sec).
    /// * `etp1` - actual transpiration rate (mm / sec).
    ///
    /// The following global variables are referenced here:
    ///
    /// AirTemp, albedo, CloudCoverRatio, CloudTypeCorr, FieldCapacity,
    /// MulchTemp, MulchTranSW, PlantHeight, Radiation, RelativeHumidity,
    /// rracol, SitePar, thad, VolWaterContent, WindSpeed.
    ///
    /// The following global variables are set here:
    /// bEnd, SoilTemp, FoliageTemp, MulchTemp.
    unsafe fn soil_energy_balance(
        &mut self,
        ihr: usize,
        k: usize,
        bMulchon: bool,
        ess: f64,
        etp1: f64,
    ) -> Result<(), Cotton2KError> {
        // Stefan-Boltsman constant.
        const stefa1: f64 = 1.38e-12;
        // Ratio of wind speed under partial canopy cover.
        const wndfac: f64 = 0.60;
        // proportion of short wave radiation (on fully shaded soil surface)intercepted by the canopy.
        const cswint: f64 = 0.75;
        // Set initial values
        // fraction of shaded soil area
        let sf = 1. - self.rracol[k];
        // air temperature, K
        let thet = AirTemp[ihr] + 273.161;
        // soil surface temperature, K
        let mut so = SoilTemp[0][k];
        // 2nd soil layer temperature, K
        let mut so2 = SoilTemp[1][k];
        // 3rd soil layer temperature, K
        let mut so3 = SoilTemp[2][k];
        // Compute soil surface albedo (based on Horton and Chung, 1991):
        // albedo of the soil surface
        let ag = if VolWaterContent[0][k] <= thad[0] {
            SitePar[15]
        } else if VolWaterContent[0][k] >= FieldCapacity[0] {
            SitePar[16]
        } else {
            SitePar[16]
                + (SitePar[15] - SitePar[16]) * (FieldCapacity[0] - VolWaterContent[0][k])
                    / (FieldCapacity[0] - thad[0])
        };
        //  ****   SHORT WAVE RADIATION ENERGY BALANCE   ****
        // Division by 41880 (= 698 * 60) converts from Joules per sq m to
        // langley (= calories per sq cm) Or: from Watt per sq m to langley per sec.
        // Modify incoming short wave radiation to mulched soil surface.
        let rzero = Radiation[ihr] / 41880.; // short wave (global) radiation (ly / sec).
        let rss0 = rzero * (1. - sf * cswint); // global radiation after passing through canopy
        let mut tm; // temperature of mulch (K)
        let rsup; // global radiation reflected up to the vegetation
        let mut rsm = 0.; // global radiation absorbed by mulch
        let rss; // global radiation absorbed by soil surface
        if bMulchon {
            tm = MulchTemp[k];
            // Assume all non transfered radiation is absorbed by mulch
            rss = rss0 * MulchTranSW * (1. - ag); // absorbed by soil surface
            rsm = rss0 * (1. - MulchTranSW) // absorbed by mulch
                      + rss0 * MulchTranSW * ag * (1. - MulchTranSW) // reflected up from soil surface and absorbed by mulch
                      ;
            rsup = rss0 * MulchTranSW * ag * MulchTranSW; // reflected up from soil surface through mulch
        } else {
            tm = 0.;
            rss = rss0 * (1. - ag); // absorbed by soil surface
            rsup = rss0 * ag; // reflected up from soil surface
        }
        //   ****   LONG WAVE RADIATION EMITTED FROM SKY    ****
        // air vapor pressure, KPa.
        let vp = 0.01 * RelativeHumidity[ihr] * vapor_pressure(AirTemp[ihr]);
        // sky emissivity from clear portions of the sky.
        let ea0 = clearskyemiss(vp, thet);
        // incoming long wave radiation (ly / sec).
        let rlzero = (ea0 * (1. - CloudCoverRatio[ihr]) + CloudCoverRatio[ihr]) * stefa1 * thet.powi(4)
                - CloudTypeCorr[ihr] / 41880. // CloudTypeCorr converted from W m-2 to ly sec-1.
                ;
        // Set initial values of canopy temperature and air temperature in canopy.
        let mut tv = 0.; // temperature of plant foliage (K)
        let mut tafk; // temperature (K) of air inside the canopy.
        if sf < 0.05 {
            // no vegetation
            tv = thet;
        }
        // Wind velocity in canopy is converted to cm / s.
        let wndhr = WindSpeed[ihr] * 100.;
        let mut rocp; // air density * specific heat at constant pressure = 0.24 * 2 * 1013 / 5740 divided by tafk.
        let mut c2 = 0.; // multiplier for sensible heat transfer (at plant surface).
        let mut rsv = 0.; // global radiation absorbed by the vegetation
        if sf >= 0.05 {
            // a shaded soil column
            // vegetation temperature
            tv = FoliageTemp[k];
            // Short wave radiation intercepted by the canopy:
            rsv = rzero * (1. - albedo[ihr]) * sf * cswint  //  from above
                      + rsup * (1. - albedo[ihr]) * sf * cswint //  reflected from soil surface
                      ;
            // Air temperature inside canopy is the average of soil, air, and plant temperatures,
            // weighted by 0.1, 0.3, and 0.6, respectively. In case of mulch, mulch temperature replaces soil temperature.
            tafk = (1. - sf) * thet
                + sf * (0.1 * if bMulchon { tm } else { so } + 0.3 * thet + 0.6 * tv);
            // Call SensibleHeatTransfer() to compute sensible heat transfer coefficient. Factor 2.2 for sensible heat transfer: 2 sides of leaf plus stems and petioles.
            // sensible heat transfer coefficient for soil
            let varcc = SensibleHeatTransfer(tv, tafk, PlantHeight, wndhr); // canopy to air
            if bEnd {
                return;
            }
            rocp = 0.08471 / tafk;
            c2 = 2.2 * sf * rocp * varcc;
        }
        // counter of iteration number
        let mut menit: usize = 0;
        // previous value of soil surface temperature
        let mut soold;
        // previous value of vegetation temperature
        let mut tvold = tv;
        // previous value of temperature of mulch (K)
        let mut tmold = tm;
        // Starting iterations for soil, mulch and canopy energy balance
        loop {
            soold = so;
            let wndcanp = (1. - sf * (1. - wndfac)) * wndhr; // estimated wind speed under canopy
            if bMulchon {
                // This section executed for mulched columns: call SoilMulchBalance() for soil / mulch interface.
                tmold = tm;
                bEnd = SoilMulchBalance(
                    ihr as i32, k as i32, rlzero, rsm, rss, sf, &mut so, &mut so2, &mut so3, thet,
                    &mut tm, tv, wndcanp,
                )?;
                if bEnd {
                    return;
                }
            } else {
                // This section executed for non-mulched columns
                // Call SensibleHeatTransfer() to compute sensible heat transfer for soil surface to air
                tafk = (1. - sf) * thet + sf * (0.1 * so + 0.3 * thet + 0.6 * tv);
                // sensible heat transfer coefficientS for soil
                let varc = SensibleHeatTransfer(so, tafk, 0., wndcanp);
                if bEnd {
                    return;
                }
                rocp = 0.08471 / tafk;
                // multiplier for computing sensible heat transfer soil to air.
                let hsg = rocp * varc;
                // Call SoilSurfaceBalance() for energy balance in soil surface/air interface.
                SoilSurfaceBalance(
                    ihr as i32, k as i32, ess, rlzero, rss, sf, hsg, &mut so, &mut so2, &mut so3,
                    thet, 0., tv,
                );
                if bEnd {
                    return;
                }
            }
            if sf >= 0.05 {
                // This section executed for shaded columns only.
                tvold = tv;
                // Compute canopy energy balance for shaded columns
                CanopyBalance(
                    ihr as i32, k as i32, etp1, rlzero, rsv, c2, sf, so, thet, tm, &mut tv,
                );
                // Increment the number of iterations.
                menit += 1;
                if menit > 10 {
                    // The following is used to reduce fluctuations.
                    so = 0.5 * (so + soold);
                    if bMulchon {
                        tm = 0.5 * (tm + tmold);
                    }
                    tv = 0.5 * (tv + tvold);
                }
                if menit > 30 {
                    // If more than 30 iterations are needed - stop simulation.
                    bEnd = true;
                    return;
                }
            }
            if !((tv - tvold).abs() > 0.05
                || (so - soold).abs() > 0.05
                || (tm - tmold).abs() > 0.05)
            {
                break;
            }
        }
        // After convergence - set global variables for the following temperatures:
        if sf >= 0.05 {
            FoliageTemp[k] = tv;
        }
        SoilTemp[0][k] = so;
        SoilTemp[1][k] = so2;
        SoilTemp[2][k] = so3;
        MulchTemp[k] = if bMulchon { tm } else { 0. };
        Ok(())
    }

    /// This function computes heat flux in one direction between soil cells.
    /// It is called from SoilTemperature(), and calls ThermalCondSoil() and
    /// [Soil::heat_balance()].
    ///
    /// Note: the units are: thermal conductivity = cal cm-1 s-1 oC-1;
    ///                      heat capacity = cal cm-3 oC-1;
    ///                      thermal diffusivity = cm2 s-1;
    ///                      ckx and cky are dimensionless;
    ///
    /// The following global variables are referenced here:
    ///
    /// dl, HeatCapacitySoilSolid, maxl, PoreSpace, VolWaterContent, wk.
    ///
    /// The following global variable is set here:    SoilTemp.
    ///
    /// The following input arguments are used in this function:
    ///
    /// dlt - time (seconds) of one iteration.
    /// iv -  = 1 for vertical flux, = 0 for horizontal flux.
    /// layer - soil layer number.
    /// n0 - number of layer or column of this array
    /// nn - number of soil cells in the array.
    unsafe fn heat_flux(&mut self, dlt: f64, iv: i32, nn: usize, layer: usize, n0: usize) {
        // Constant parameters:
        // weighting factor for the implicit method of computation.
        const beta1: f64 = 0.90;
        // heat capacity of air (cal cm-3 oC-1).
        const ca: f64 = 0.0003;
        // Set soil layer number l (needed to define HeatCapacitySoilSolid, PoreSpace, ThermalCondSoil).
        // Compute for each soil cell the heat capacity and heat diffusivity.
        let mut l = layer; // soil layer number.
        let mut q1: [f64; 40] = [0.; 40]; // array of water content.
        let mut asoi: [f64; 40] = [0.; 40]; // array of thermal diffusivity of soil cells (cm2 s-1).
        for i in 0..nn {
            if iv == 1 {
                l = i;
                q1[i] = VolWaterContent[i][n0];
                self.ts1[i] = SoilTemp[i][n0];
                self.dz[i] = dl[i];
            } else {
                q1[i] = VolWaterContent[n0][i];
                self.ts1[i] = SoilTemp[n0][i];
                self.dz[i] = wk[i];
            }
            self.hcap[i] = HeatCapacitySoilSolid[l] + q1[i] + (PoreSpace[l] - q1[i]) * ca;
            asoi[i] = ThermalCondSoil(q1[i], self.ts1[i], l as i32) / self.hcap[i];
        }
        // The numerical solution of the flow equation is a combination of the implicit method (weighted by beta1) and the explicit method (weighted by 1-beta1).
        let mut dltt; // computed time step required.
        let mut avdif: [f64; maxl] = [0.; maxl]; // average thermal diffusivity between adjacent cells.
        let mut dy: [f64; maxl] = [0.; maxl]; // array of distances between centers of adjacent cells (cm).
        let mut dltmin = dlt; // minimum time step for the explicit solution.
        for i in 1..nn {
            // Compute average diffusivities avdif between layer i and the previous (i-1),
            //  and dy(i), distance (cm) between centers of layer i and the previous
            //  (i-1)
            avdif[i] = (asoi[i] + asoi[i - 1]) / 2.;
            dy[i] = (self.dz[i - 1] + self.dz[i]) / 2.;
            //     Determine the minimum time step required for the explicit
            //     solution.
            dltt = 0.2 * dy[i] * self.dz[i] / avdif[i] / (1. - beta1);
            if dltt < dltmin {
                dltmin = dltt;
            }
        }
        //     Use time step of dlt1 seconds, for iterx iterations
        let mut iterx = (dlt / dltmin) as usize; // computed number of iterations.
        if dltmin < dlt {
            iterx += 1;
        }
        let dlt1 = dlt / iterx as f64; // computed time (seconds) of an iteration.
                                       // start iterations. Store temperature data in array ts0. count iterations.
        for _ in 0..iterx {
            for i in 0..nn {
                self.ts0[i] = self.ts1[i];
                if iv == 1 {
                    l = i;
                }
                asoi[i] = ThermalCondSoil(q1[i], self.ts1[i], l as i32) / self.hcap[i];
                if i > 0 {
                    avdif[i] = (asoi[i] + asoi[i - 1]) / 2.;
                }
            }
            self.numiter += 1;
            // The solution of the simultaneous equations in the implicit method alternates between
            // the two directions along the arrays. The reason for this is because
            // the direction of the solution may cause some cumulative bias. The
            // counter numiter determines the direction of the solution.
            let mut cau: [f64; maxl] = [0.; maxl];
            let mut dau: [f64; maxl] = [0.; maxl]; // arrays used for the implicit numerical solution.
            let mut ckx = 0.;
            let mut cky = 0.; // nondimensional diffusivities to next and previous layers.
            let mut vara;
            let mut varb; // used for computing the implicit solution.
            if (self.numiter % 2) == 0 {
                // 1st direction of computation, for an even iteration number:
                dau[0] = 0.;
                cau[0] = self.ts1[0];
                // Loop from the second to the last but one soil cells. Compute nondimensional diffusivities to next and previous layers.
                for i in 1..(nn - 1) {
                    ckx = avdif[i + 1] * dlt1 / (self.dz[i] * dy[i + 1]);
                    cky = avdif[i] * dlt1 / (self.dz[i] * dy[i]);
                    // Correct value of layer 1 for explicit heat movement to/from layer 2
                    if i == 1 {
                        cau[0] = self.ts1[0]
                            - (1. - beta1) * (self.ts1[0] - self.ts1[1]) * cky * self.dz[1]
                                / self.dz[0];
                    }
                    vara = 1. + beta1 * (ckx + cky) - beta1 * ckx * dau[i - 1];
                    dau[i] = beta1 * cky / vara;
                    varb = self.ts1[i]
                        + (1. - beta1)
                            * (cky * self.ts1[i - 1] + ckx * self.ts1[i + 1]
                                - (cky + ckx) * self.ts1[i]);
                    cau[i] = (varb + beta1 * ckx * cau[i - 1]) / vara;
                }
                // Correct value of last layer (nn-1) for explicit heat movement to/from layer nn-2
                self.ts1[nn - 1] = self.ts1[nn - 1]
                    - (1. - beta1) * (self.ts1[nn - 1] - self.ts1[nn - 2]) * ckx * self.dz[nn - 2]
                        / self.dz[nn - 1];
                // Continue with the implicit solution
                for i in (0..nn - 1).rev() {
                    self.ts1[i] = dau[i] * self.ts1[i + 1] + cau[i];
                }
            } else {
                // Alternate direction of computation for odd iteration number
                dau[nn - 1] = 0.;
                cau[nn - 1] = self.ts1[nn - 1];
                for i in (0..nn - 1).rev() {
                    ckx = avdif[i + 1] * dlt1 / (self.dz[i] * dy[i + 1]);
                    cky = avdif[i] * dlt1 / (self.dz[i] * dy[i]);
                    if i == nn - 2 {
                        cau[nn - 1] = self.ts1[nn - 1]
                            - (1. - beta1)
                                * (self.ts1[nn - 1] - self.ts1[nn - 2])
                                * ckx
                                * self.dz[nn - 2]
                                / self.dz[nn - 1];
                    }
                    vara = 1. + beta1 * (ckx + cky) - beta1 * cky * dau[i + 1];
                    dau[i] = beta1 * ckx / vara;
                    varb = self.ts1[i]
                        + (1. - beta1)
                            * (ckx * self.ts1[i + 1] + cky * self.ts1[i - 1]
                                - (cky + ckx) * self.ts1[i]);
                    cau[i] = (varb + beta1 * cky * cau[i + 1]) / vara;
                }
                self.ts1[0] = self.ts1[0]
                    - (1. - beta1) * (self.ts1[0] - self.ts1[1]) * cky * self.dz[1] / self.dz[0];
                for i in 1..nn {
                    self.ts1[i] = dau[i] * self.ts1[i - 1] + cau[i];
                }
            }
            // Call self.heat_balance to correct quantitative deviations caused by the imlicit part of the solution.
            self.heat_balance(nn);
        }
        // Set values of SoiTemp
        for i in 0..nn {
            if iv == 1 {
                SoilTemp[i][n0] = self.ts1[i];
            } else {
                SoilTemp[n0][i] = self.ts1[i];
            }
        }
    }
    /// Checks and corrects the heat balance in the soil soil cells, within
    /// a soil layer. It is called by function SoilHeatFlux() only for
    /// horizontal flux.
    ///
    /// The implicit part of the solution may cause some deviation in the total
    /// heat sum to occur. This module corrects the heat balance if the sum of
    /// absolute deviations is not zero, so that the total amount of heat in
    /// the array does not change. The correction is proportional to the
    /// difference between the previous and present heat amounts.
    ///
    /// The following arguments are referenced here:
    ///
    /// nn - the number of soil cells in this layer or column.
    ///
    /// The following global or file scope variables are referenced:
    ///
    /// dz, hcap, ts0
    ///
    /// The following file scope variable is set:    ts1
    fn heat_balance(&mut self, nn: usize) {
        let mut dabs = 0.; // Sum of absolute value of differences in heat content in
                           // the array between beginning and end of this time step.
        let mut dev = 0.; // Sum of differences of heat amount in soil.
        for i in 0..nn {
            dev += self.dz[i] * self.hcap[i] * (self.ts1[i] - self.ts0[i]);
            dabs += (self.ts1[i] - self.ts0[i]).abs();
        }
        if dabs > 0. {
            for i in 0..nn {
                self.ts1[i] -=
                    (self.ts1[i] - self.ts0[i]).abs() * dev / (dabs * self.dz[i] * self.hcap[i]);
            }
        }
    }
}

/*
                           References.

     Benjamin, J.G., Ghaffarzadeh, M.R. and Cruse, R.M. 1990.
  Coupled water and heat transport in ridged soils. Soil Sci. Soc.
  Am. J. 54:963-969.

     Chen, J. 1984. Uncoupled multi-layer model for the transfer of
  sensible and latent heat flux densities from vegetation. Boundary-
  Layer Meteorology 28:213-225.

     Chen, J. 1985. A graphical extrapolation method to determine
  canopy resistance from measured temperature and humidity profiles
  above a crop canopy. Agric. For. Meteorol. 37:75-88.

     Clothier, B.E., Clawson, K.L., Pinter, P.J.Jr., Moran, M.S.,
  Reginato, R.J. and Jackson, R.D. 1986. Estimation of soil heat flux
  from net radiation during the growth of alfalfa. Agric. For.
  Meteorol. 37:319-329.

     Costello, T.A. and Braud, H.J. Jr. 1989. Thermal diffusivity
  of soil by nonlinear regression analysis of soil temperature data.
  Trans. ASAE 32:1281-1286.

     De Vries, D.A. 1963. Thermal properties of soils. In: W.R. Van
  Wijk (ed) Physics of plant environment, North Holland, Amsterdam,
  pp 210-235.

     Deardorff, J.W. 1978. Efficient prediction of ground surface
  temperature and moisture with inclusion of a layer of vegetation.
  J. Geophys. Res. 83 (C4):1889-1903.

     Dong, A., Prashar, C.K. and Grattan, S.R. 1988. Estimation of
  daily and hourly net radiation. CIMIS Final Report June 1988, pp.
  58-79.

     Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling
  diurnal patterns of air temperature, radiation, wind speed and
  relative humidity by equations from daily characteristics.
  Agricultural Systems 51:377-393.

     Hadas, A. 1974. Problem involved in measuring the soil thermal
  conductivity and diffusivity in a moist soil. Agric. Meteorol.
  13:105-113.

     Hadas, A. 1977. Evaluation of theoretically predicted thermal
  conductivities of soils under field and laboratory conditions. Soil
  Sci. Soc. Am. J. 41:460-466.

     Hanks, R.J., Austin, D.D. and Ondrechen, W.T. 1971. Soil
  temperature estimation by a numerical method. Soil Sci. Soc. Am.
  Proc. 35:665-667.

     Hares, M.A. and Novak, M.D. 1992. Simulation of surface energy
  balance and soil temperature under strip tillage: I. Model
  description. Soil Sci. Soc. Am. J. 56:22-29.

     Hares, M.A. and Novak, M.D. 1992. Simulation of surface energy
  balance and soil temperature under strip tillage: II. Field test.
  Soil Sci. Soc. Am. J. 56:29-36.

     Horton, E. and Wierenga, P.J. 1983. Estimating the soil heat
  flux from observations of soil temperature near the surface. Soil
  Sci. Soc. Am. J. 47:14-20.

     Horton, E., Wierenga, P.J. and Nielsen, D.R. 1983. Evaluation
  of methods for determining apparent thermal diffusivity of soil
  near the surface. Soil Sci. Soc. Am. J. 47:25-32.

     Horton, R. 1989. Canopy shading effects on soil heat and water
  flow. Soil Sci. Soc. Am. J. 53:669-679.

     Horton, R., and Chung, S-O, 1991. Soil Heat Flow. Ch. 17 in: Hanks,
  J., and Ritchie, J.T., (Eds.) Modeling Plant and Soil Systems. Am. Soc.
  Agron., Madison, WI, pp 397-438.

     Iqbal, M. 1983. An Introduction to Solar Radiation. Academic
  Press.

     Kimball, B.A., Jackson, R.D., Reginato, R.J., Nakayama, F.S.
  and Idso, S.B. 1976. Comparison of field-measured and calculated
  soil heat fluxes. Soil Sci. Soc. Am. J. 40:18-28.

     Lettau, B. 1971. Determination of the thermal diffusivity in
  the upper layers of a natural ground cover. Soil Sci. 112:173-177.

     Mahrer, Y. 1979. Prediction of soil temperatures of a soil
  mulched with transparent polyethylene. J. Appl. Meteorol. 18:1263-
  1267.

     Mahrer, Y. 1980. A numerical model for calculating the soil
  temperature regime under transparent polyethylene mulches. Agric.
  Meteorol. 22:227-234.

     Mahrer, Y., Naot, O., Rawitz, E. and Katan, J. 1984.
  Temperature and moisture regimes in soils mulched with transparent
  polyethylene. Soil Sci. Soc. Amer. J. 48:362-367.

     Monin, A.S. 1973. Boundary layers in planetary atmospheres. In:
  P. Morrel (ed.), Dynamic meteorology, D. Reidel Publishing Company,
  Boston, pp. 419-458.

     Spitters, C.J.T., Toussaint, H.A.J.M. and Goudriaan, J. 1986.
  Separating the diffuse and direct component of global radiation and
  its implications for modeling canopy photosynthesis. Part I.
  Components of incoming radiation. Agric. For. Meteorol. 38:217-229.

     Wierenga, P.J. and de Wit, C.T. 1970. Simulation of heat flow
  in soils. Soil Sci. Soc. Am. Proc. 34:845-848.

     Wierenga, P.J., Hagan, R.M. and Nielsen, D.R. 1970. Soil
  temperature profiles during infiltration and redistribution of cool
  and warm irrigation water. Water Resour. Res. 6:230-238.

     Wierenga, P.J., Nielsen, D.R. and Hagan, R.M. 1969. Thermal
  properties of soil based upon field and laboratory measurements.
  Soil Sci. Soc. Am. Proc. 33:354-360.
*/

/// This function solves the energy balance equations at the interface of
/// the soil surface and the plastic mulch cover and computes the resulting
/// temperatures of the soil surface and of the plastic mulch.
/// It is called from [EnergyBalance()], on each time step and for each
/// soil column, if this column is covered with a plastic mulch.  It calls functions
//  SensibleHeatTransfer(), SoilSurfaceBalance() and MulchSurfaceBalance().
///
/// Units for all energy fluxes are: cal cm-2 sec-1.
///
/// If the return value is true, it means there was an error and simulation will end.
///
/// The following global variables are referenced here:
///
/// bEnd, Daynum, MulchTranLW, .
///
/// The following arguments are set in this function:
///
/// so - temperature of soil surface.
/// so2 - temperature of soil 2nd layer
/// so3 - temperature of soil 3rd layer
/// tm - temperature of plastic mulch (K). When tm = 0 there is no plastic mulch.
///
/// The following arguments are referenced in this function:
///
/// ihr - the time in hours.
/// k - soil column number.
/// rlzero - incoming long wave radiation (ly / sec).
/// rsm - global radiation absorbed by mulch
/// rss - global radiation absorbed by soil surface
/// sf - fraction of shaded soil area
/// thet - air temperature (K).
/// tv - temperature of plant canopy (K).
/// wndcanp - estimated wind speed under canopy
fn SoilMulchBalance(
    ihr: i32,
    k: i32,
    rlzero: f64,
    rsm: f64,
    rss: f64,
    sf: f64,
    so: &mut f64,
    so2: &mut f64,
    so3: &mut f64,
    thet: f64,
    tm: &mut f64,
    tv: f64,
    wndcanp: f64,
) -> Result<bool, Cotton2KError> {
    // Constant variables:
    // emissivity of the foliage surface.
    const ef: f64 = 0.95;
    // emissivity of the soil surface.
    const eg: f64 = 0.95;
    // stefan-boltsman constant.
    const stefa1: f64 = 1.38e-12;
    // Compute long wave radiation reaching the surface mulch from above, and the air temperature above it.
    let rlsp0; // long wave radiation reaching the mulch from above.
    let tafk; //  temperature (K) of air inside the canopy.
    if sf > 0.05
    // shaded column
    {
        rlsp0 = (1. - sf) * (1. - unsafe {MulchTranLW}) * rlzero  // from sky in unshaded segment
+ sf * (1. - unsafe {MulchTranLW}) * ef * stefa1 * tv.powi(4); // from foliage in shaded segment
        tafk = (1. - sf) * thet + sf * (0.1 * *tm + 0.3 * thet + 0.6 * tv);
    } else
    // unshaded column
    {
        rlsp0 = (1. - unsafe { MulchTranLW }) * rlzero;
        tafk = thet;
    }
    // rls5 is the multiplier of tm**4 for emitted long wave radiation from mulch, sum of upward and downward emittance.
    // multiplier for emitted long wave radiation from mulch,
    let rls5 = 2. * (1. - unsafe { MulchTranLW }) * stefa1;
    // Call SensibleHeatTransfer() to compute sensible heat transfer between plastic mulch and air sensible heat transfer coefficients for mulch to air (before multiplying by ROCP).
    let varcm = unsafe { SensibleHeatTransfer(*tm, tafk, 0., wndcanp) };
    if unsafe { bEnd } {
        return Ok(true);
    }
    // multiplier for computing sensible heat transfer from soil to mulch.
    let mut hsgm;
    // air density * specific heat at constant pressure [= 0.24 * 2 * 1013 / (5740 * tk) ]
    let mut rocp = 0.08471 / tafk;
    // multiplier for computing sensible heat transfer from mulch to air.
    let hsgp = rocp * varcm;
    // Compute sensible heat transfer between plastic mulch and soil surface
    if *tm > 0. {
        rocp = 0.08471 / *tm;
    }
    let mut mtnit = 0; // counter for numbet of iterations
    let mut soold1; // previous value of temperature of soil surface (k)
    let mut tmold1; // previous value of temperature of mulch (k)

    loop {
        soold1 = *so;
        tmold1 = *tm;
        // Energy balance for soil surface (mulch interface)
        hsgm = 2. * rocp * (*so - *tm).abs();
        unsafe {
            SoilSurfaceBalance(
                ihr, k, 0., rlzero, rss, sf, hsgm, so, so2, so3, thet, *tm, tv,
            )
        };
        if unsafe { bEnd } {
            return Ok(true);
        }
        // Add Long wave radiation reaching the mulch from the soil: total long wave radiation reaching the mulch.
        let rlsp = rlsp0 + (1. - unsafe { MulchTranLW }) * eg * stefa1 * so.powi(4);
        // Energy balance for mulch (soil and air interface)
        hsgm = 2. * rocp * (*so - *tm).abs();
        unsafe { MulchSurfaceBalance(ihr, k, rlsp, rls5, rsm, sf, hsgp, hsgm, *so, thet, tm, tv) };
        if unsafe { bEnd } {
            return Ok(true);
        }
        // Check number of iterations - do not exceed 30 iterations.
        mtnit += 1;
        if mtnit > 8 {
            *so = (*so + soold1) / 2.;
            *tm = (*tm + tmold1) / 2.;
        }
        if mtnit > 30 {
            return Err(Cotton2KError {
                level: 0,
                message: String::from("Infinite loop in SoilMulchBalance(). Abnormal stop!! \n"),
            });
        }
        if (*tm - tmold1).abs() <= 0.05 && (*so - soold1).abs() <= 0.05 {
            return Ok(false);
        }
    }
}
