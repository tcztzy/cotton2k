use crate::{
    bEnd, dl, es1hour, es2hour, isw, nk, nl, rracol, thad, wk, ActualSoilEvaporation,
    ActualTranspiration, Clim, CumEvaporation, DayEndMulch, DayOfSimulation, DayPlant, DayStart,
    DayStartMulch, Daynum, DeepSoilTemperature, EnergyBalance, FoliageTemp, Kday, MulchIndicator,
    MulchTemp, PlantRowColumn, PredictEmergence, ReferenceETP, ReferenceTransp, RowSpace,
    Scratch21, SitePar, SoilHeatFlux, SoilTemp, SoilTempDailyAvrg, VolWaterContent,
};

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
pub unsafe fn SoilTemperature() {
    if Daynum <= DayStart {
        SoilTemperatureInit();
    }
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
        kk = nk;
        for k in 0..nk as usize {
            if rracol[k] <= 0.99 {
                shadetot += 1. - rracol[k];
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
        let mut tmav = 0.; // average mulch temperature.
        let mut kmulch = 0u32; // number of soil columns covered with mulch.

        // Compute vertical transport for each column
        for k in 0..kk as usize {
            //     Set SoilTemp for the lowest soil layer.
            SoilTemp[(nl - 1) as usize][k] = DeepSoilTemperature;
            //     Compute transpiration from each column, weighted by its
            //     relative shading.

            // actual hourly transpiration (mm s-1) for a column.
            let etp1 = if shadeav > 0.000001 {
                etp0 * (1. - rracol[k]) / shadeav
            } else {
                0.
            };
            //     Check if mulch is on for this date and for this column.
            // is true if this column is covered with plastic mulch now, false if not.
            let bMulchon = if MulchIndicator == 0 || Daynum < DayStartMulch || Daynum > DayEndMulch
            {
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
                    if k >= (PlantRowColumn - 1) as usize && k <= (PlantRowColumn + 2) as usize {
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
                let mut escol1k = es1hour[ihr] * rracol[k] + es2hour[ihr];
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
            EnergyBalance(ihr as i32, k as i32, bMulchon, ess, etp1);
            if bEnd {
                return;
            }
            if bMulchon {
                tmav += MulchTemp[k] - 273.161;
                kmulch += 1;
            }
        }

        if kmulch > 0 {
            tmav /= kmulch as f64;
        }
        // Compute soil temperature flux in the vertical direction.
        // Assign iv = 1, layer = 0, nn = nl.
        let mut iv = 1; // indicates vertical (=1) or horizontal (=0) flux.
        let mut nn = nl; // number of array members for heat flux.
        let mut layer = 0; // soil layer number
        let mut tsolav: [f64; 40] = [0.; 40]; // hourly average soil temperature C, of a soil layer.
        for k in 0..kk {
            // Loop over kk columns, and call SoilHeatFlux().
            SoilHeatFlux(dlt, iv, nn, layer, k);
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
            nn = nk;
            for l in 0..nl {
                layer = l;
                SoilHeatFlux(dlt, iv, nn, layer, l);
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
        // Compute average temperature of foliage, in degrees C. The average is weighted by the canopy shading of each column, only columns which are shaded 5% or more by canopy are used.
        let mut tfc = 0.; // average foliage temperature, weighted by shading in each column
        let mut shading = 0.; // sum of shaded area in all shaded columns, used to compute TFC
        for k in 0..nk as usize {
            if rracol[k] <= 0.95 {
                tfc += (FoliageTemp[k] - 273.161) * (1. - rracol[k]);
                shading += 1. - rracol[k];
            }
        }
        if shading >= 0.01 {
            tfc /= shading;
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

//     This function is called from SoilTemperature() at the start of the
//     simulation. It sets
//  initial values to soil and canopy temperatures.
//     The following arguments are set in this function:
//  jt1, jt2 - input of start and stop of output of soil temperatures.
//
//     The following global variables are referenced here:
//  Clim (structure), DayFinish, Daynum, DayStart, nl,
//  SitePar.
//     The following global variables are set here:
//  DeepSoilTemperature, SoilTemp.
pub unsafe fn SoilTemperatureInit() {
    //     If there is an output flag for soil temperatures, an error message
    //     pops up for defining
    //  the start and stop dates for this output.
    //     Compute initial values of soil temperature: It is assumed that at the
    //     start of simulation
    //  the temperature of the first soil layer (upper boundary) is equal to the
    //  average air temperature of the previous five days (if climate data not
    //  available - start from first climate data).
    //     NOTE: For a good simulation of soil temperature, it is recommended to
    //     start simulation at
    //  least 10 days before planting date. This means that climate data should
    //  be available for this period. This is especially important if emergence
    //  date has to be simulated.
    let mut idd = Daynum - 4 - DayStart; // number of days minus 4 from start of simulation.
    if idd < 0 {
        idd = 0;
    }
    let mut tsi1 = 0.; // Upper boundary (surface layer) initial soil temperature, C.
    for i in idd as usize..(idd + 5) as usize {
        tsi1 += Clim[i].Tmax + Clim[i].Tmin;
    }
    tsi1 = tsi1 / 10.;
    //     The temperature of the last soil layer (lower boundary) is computed
    //     as a sinusoidal function
    //  of day of year, with site-specific parameters.
    DeepSoilTemperature = SitePar[9]
        + SitePar[10] * (2. * std::f64::consts::PI * (Daynum as f64 - SitePar[11]) / 365.).sin();
    //     SoilTemp is assigned to all columns, converted to degrees K.
    tsi1 += 273.161;
    DeepSoilTemperature += 273.161;
    for l in 0..nl as usize {
        //     The temperatures of the other soil layers are linearly
        //     interpolated.
        //  tsi = computed initial soil temperature, C, for each layer
        let tsi = ((nl as usize - l - 1) as f64 * tsi1 + l as f64 * DeepSoilTemperature)
            / (nl - 1) as f64;
        for k in 0..nk as usize {
            SoilTemp[l][k] = tsi;
        }
    }
}
