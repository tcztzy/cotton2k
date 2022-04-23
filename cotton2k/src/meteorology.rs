use crate::{
    albedo, bPollinSwitch, es1hour, es2hour, fmin, iyear, AirTemp, AvrgDailyTemp,
    ClayVolumeFraction, Clim, CloudCoverRatio, CloudTypeCorr, DayLength, DayOfSimulation, DayStart,
    DayTimeTemp, Daynum, DewPointTemp, Elevation, GetFromClim, Irrig, LastDayWeatherData, Latitude,
    LeapYear, Longitude, NightTimeTemp, NumIrrigations, Radiation, ReferenceETP, ReferenceTransp,
    RelativeHumidity, Rn, SandVolumeFraction, Scratch21, SitePar, WindSpeed, CLIMATE_METRIC_IRRD,
    CLIMATE_METRIC_RAIN, CLIMATE_METRIC_TDEW, CLIMATE_METRIC_TMAX, CLIMATE_METRIC_TMIN,
    CLIMATE_METRIC_WIND,
};
/// daily declination angle, in radians.
static mut declination: f64 = 0.;
/// time of solar noon, hours.
static mut SolarNoon: f64 = 0.;
/// time of sunrise, hours.
static mut sunr: f64 = 0.;
/// time of sunset, hours.
static mut suns: f64 = 0.;
/// extraterrestrial radiation, W / m2.
static mut tmpisr: f64 = 0.;
//     All daily weather data, read from input file, are stored in the
//     structure:
//         Clim[400];    --  defined in global
//     The values are extracted from this structure by function GetFromClim(),
//     see
//  file "GeneralFunctions.cpp"

//     The function DayClim() is called daily from SimulateThisDay(). It calls
//     the
//  the following functions:
//     ComputeDayLength(), GetFromClim(), SimulateRunoff(),
//     AverageAirTemperatures(), dayrad(), daytemp(), EvapoTranspiration(),
//     tdewhour(), dayrh(), daywnd()
//         Global variables referenced:
//    Daynum, DayFinish, DayStart, declination, LastDayWeatherData, Latitude,
//    pi, SitePar
//         Global variables set:
//    AirTemp, bPollinSwitch, DewPointTemp, Radiation, RelativeHumidity,
//    WindSpeed
pub unsafe fn DayClim() {
    //     Compute day length and related variables:
    ComputeDayLength();
    //
    let xlat = Latitude * std::f64::consts::PI / 180.; // latitude converted to radians.
    let cd = xlat.cos() * declination.cos(); // amplitude of the sine of the solar height.
    let sd = xlat.sin() * declination.sin(); // seasonal offset of the sine of the solar height.
                                             // The computation of the daily integral of global radiation (from sunrise to sunset) is based on Spitters et al. (1986).
    const c11: f64 = 0.4; // constant parameter.
    let radsum =            // daily radiation integral.
    if (sd / cd).abs() >= 1. {
        //  arctic circle
        0.
    } else {
        //     dsbe is the integral of sinb * (1 + c11 * sinb) from sunrise to
        //     sunset,
        let dsbe = (-sd / cd).acos() * 24. / std::f64::consts::PI * (sd + c11 * sd * sd + 0.5 * c11 * cd * cd) + 12. * (cd * (2. + 3. * c11 * sd)) * (1. - (sd / cd) * (sd / cd)).sqrt() / std::f64::consts::PI;
        //     The daily radiation integral is computed for later use in
        //     function Radiation.
        //  Daily radiation intedral is converted from langleys to Watt m-2, and
        //  divided by dsbe.
        //      11.630287 = 1000000 / 3600 / 23.884
        GetFromClim(CLIMATE_METRIC_IRRD, Daynum) * 11.630287 / dsbe
    };
    // Set 'pollination switch' for rainy days (as in GOSSYM).
    // The amount of rain today, mm
    let mut rainToday = GetFromClim(CLIMATE_METRIC_RAIN, Daynum);
    bPollinSwitch = rainToday < 2.5;
    //     Call SimulateRunoff() only if the daily rainfall is more than 2 mm.
    //     Note: this is modified from the original GOSSYM - RRUNOFF routine. It
    //     is called here
    //  for rainfall only, but it is not activated when irrigation is applied.
    let mut runoffToday = 0.; // amount of runoff today, mm
    if rainToday >= 2. {
        runoffToday = SimulateRunoff(rainToday);
        if runoffToday < rainToday {
            rainToday -= runoffToday;
        } else {
            rainToday = 0.;
        }
        let j = Daynum - DayStart; // days from start of simulation
        Clim[j as usize].Rain = rainToday;
    }
    Scratch21[(DayOfSimulation - 1) as usize].runoff = runoffToday;
    //     Parameters for the daily wind function are now computed:
    //     Note:  SitePar[] are site specific parameters.
    let t1 = sunr + SitePar[1]; // the hour at which wind begins to blow
                                // (SitePar(1) hours after sunrise).
    let t2 = SolarNoon + SitePar[2]; // the hour at which wind speed is maximum
                                     // (SitePar(2) hours after solar noon).
    let t3 = suns + SitePar[3]; // the hour at which wind stops to blow
                                // (SitePar(3) hours after sunset).
    let wnytf = SitePar[4]; // used for estimating night time wind (from
                            // time t3 to time t1 next day).
                            //
    for ihr in 0..24 as usize {
        let ti = ihr as f64 + 0.5; // time in the middle of each hourly interval.
        let sinb = sd + cd * (std::f64::consts::PI * (ti - SolarNoon) / 12.).cos(); // sine of the solar elevation.
                                                                                    //     Compute hourly global radiation, using function dayrad.
        Radiation[ihr] = dayrad(radsum, sinb, c11);
        //     Compute hourly temperature, using function daytmp.
        AirTemp[ihr] = daytmp(ti);
        //     Compute hourly dew point temperature, using function tdewhour.
        DewPointTemp[ihr] = tdewhour(ti, AirTemp[ihr]);
        //     Compute hourly relative humidity, using function dayrh.
        RelativeHumidity[ihr] = dayrh(AirTemp[ihr], DewPointTemp[ihr]);
        //     Compute hourly wind speed, using function daywnd, and daily sum
        //     of wind.
        WindSpeed[ihr] = daywnd(
            ti,
            GetFromClim(CLIMATE_METRIC_WIND, Daynum),
            t1,
            t2,
            t3,
            wnytf,
        );
    }
    //     Compute average daily temperature, using function
    //     AverageAirTemperatures.
    AverageAirTemperatures();
    //     Compute potential evapotranspiration.
    EvapoTranspiration();
}

//     Function ComputeDayLength() computes day length, declination, time of
//  solar noon, and extra-terrestrial radiation for this day. The CIMIS
//  (California Irrigation Management Information System) algorithms are used.
//     Global variables referenced here:
//  Daynum, iyear, Latitude, Longitude, pi,
//     Global variables set here:
//  DayLength, declination
unsafe fn ComputeDayLength() {
    //     Convert day of year to corresponding angle in radians (xday). It uses
    //     function
    //  LeapYear() (see file GeneralFunctions.cpp)
    let xday = 2. * std::f64::consts::PI * (Daynum - 1) as f64 / (365 + LeapYear(iyear)) as f64;
    //     Compute declination angle for this day. The equation used here for
    //     computing it
    //  is taken from the CIMIS algorithm.
    declination = 0.006918 - 0.399912 * xday.cos() + 0.070257 * xday.sin()
        - 0.006758 * (2. * xday).cos()
        + 0.000907 * (2. * xday).sin()
        - 0.002697 * (3. * xday).cos()
        + 0.001480 * (3. * xday).sin();
    //     Compute extraterrestrial radiation in W m-2. The 'solar constant'
    //     (average
    //  value = 1367 W m-2) is corrected for this day's distance between earth
    //  and the sun. The equation used here is from the CIMIS algorithm, which
    //  is based on the work of Iqbal (1983).
    tmpisr = 1367.
        * (1.00011
            + 0.034221 * xday.cos()
            + 0.00128 * xday.sin()
            + 0.000719 * (2. * xday).cos()
            + 0.000077 * (2. * xday).sin());
    //     Time of solar noon (SolarNoon) is computed by the CIMIS algorithm,
    //  using a correction for longitude (f), and the date correction (exday).
    //     It is assumed that the time zone is "geographically correct". For
    //  example, longitude between 22.5 and 37.5 East is in time zone GMT+2,
    //  and longitude between 112.5 and 127.5 West is in time zone GMT-8.
    //     All daily times in the model are computed by this method.
    let exday = (0.000075 + 0.001868 * xday.cos()
        - 0.032077 * xday.sin()
        - 0.014615 * (2. * xday).cos()
        - 0.04089 * (2. * xday).sin())
        * 12.
        / std::f64::consts::PI;
    let st = 15. * (Longitude / 15.).floor();
    let mut f = (Longitude - st) / 15.;

    if Longitude > 0. {
        // east of Greenwich
        if f > 0.5 {
            f -= 1.;
        }
    } else {
        // west  of Greenwich
        if f < -0.5 {
            f += 1.;
        }
    }
    SolarNoon = 12. - f - exday;
    //     Compute day length, by commonly used equations, from latitude and
    //     declination of
    //  this day. Times of sunrise and of sunset are computed from solar noon
    //  and day length.
    let xlat = Latitude * std::f64::consts::PI / 180.;
    let mut ht = -xlat.tan() * declination.tan();
    if ht > 1. {
        ht = 1.; //  arctic circle
    } else if ht < -1. {
        ht = -1.;
    }
    DayLength = 2. * ht.acos() * 12. / std::f64::consts::PI;
    sunr = SolarNoon - DayLength / 2.;
    suns = sunr + DayLength;
}
//     Function dayrad() computes the hourly values of global radiation, in W
//     m-2,
//  using the measured daily total global radiation.
//     The algorithm follows the paper of Spitters et al. (1986). It assumes
//  that atmospheric transmission of radiation is lower near the margins of
//  the daylight period, because of an increase in the path length through
//  the atmosphere at lower solar heights. Radiation is therefore assumed to be
//  proportional to sinb * (1 + c11 * sinb), where the value of c11 is set as
//  0.4 .
//     Input arguments:
//        radsum - daily radiation integral.
//        sinb - sine of the solar elevation.
//        c11 - constant parameter (0.4).
//
fn dayrad(radsum: f64, sinb: f64, c11: f64) -> f64 {
    let HourlyRadiation = radsum * sinb * (1. + c11 * sinb);
    if HourlyRadiation < 0. {
        0.
    } else {
        HourlyRadiation
    }
}
//      References:
//      Spitters, C.J.T., Toussaint, H.A.J.M. and Goudriaan, J. 1986.
// Separating the diffuse and direct component of global radiation and
// its implications for modeling canopy photosynthesis. Part I.
// Components of incoming radiation. Agric. For. Meteorol. 38:217-229.
//     Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling
//  diurnal patterns of air temperature, radiation, wind speed and
//  relative humidity by equations from daily characteristics.
//  Agricultural Systems 51:377-393.

//     Function daytmp() computes and returns the hourly values of air
//     temperature,
//  using the measured daily maximum and minimum.
//     The algorithm is described in Ephrath et al. (1996). It is based on
//  the following assumptions:
//     The time of minimum daily temperature is at sunrise.
//     The time of maximum daily temperature is SitePar[8] hours after solar
//     noon. Many models assume a sinusoidal curve of the temperature during the
//     day,
//  but actual data deviate from the sinusoidal curve in the following
//  characteristic way: a faster increase right after sunrise, a near plateau
//  maximum during several hours in the middle of the day, and a rather fast
//  decrease by sunset. The physical reason for this is a more efficient mixing
//  of heated air from ground level into the atmospheric boundary layer, driven
//  by strong lapse temperature gradients buoyancy.
//     Note: ** will be used for "power" as in Fortran notation.
//     A first order approximation is
//       daytmp = tmin + (tmax-tmin) * st * tkk / (tkk + daytmp - tmin)
//  where
//       st = sin(pi * (ti - SolarNoon + dayl / 2) / (dayl + 2 * SitePar[8]))
//     Since daytmp appears on both sides of the first equation, it
//  can be solved and written explicitly as:
//      daytmp = tmin - tkk/2 + 0.5 * sqrt(tkk**2 + 4 * amp * tkk * st)
//  where the amplitude of tmin and tmax is calculated as
//      amp = (tmax - tmin) * (1 + (tmax - tmin) / tkk)
//     This ensures that temperature still passes through tmin and tmax values.
//     The value of tkk was determined by calibration as 15.
//     This algorithm is used for the period from sunrise to the time of maximum
//     temperature,
//  hmax. A similar algorithm is used for the time from hmax to sunset, but the
//  value of the minimum temperature of the next day (mint_tomorrow) is used
//  instead of mint_today.
//     Night air temperature is described by an exponentially declining curve.
//     For the time from sunset to mid-night:
//        daytmp = (mint_tomorrow - sst * exp((dayl - 24) / tcoef)
//               + (sst - mint_tomorrow) * exp((suns - ti) / tcoef))
//               / (1 - exp((dayl - 24) / tcoef))
//  where
//        tcoef is a time coefficient, determined by calibration as 4
//        sst is the sunset temperature, determined by the daytime equation as:
//        sst = mint_tomorrow - tkk / 2 + 0.5 * sqrt(tkk**2 + 4 * amp * tkk *
//        sts)
//  where
//        sts  = sin(pi * dayl / (dayl + 2 * SitePar[8]))
//        amp = (tmax - mint_tomorrow) * (1 + (tmax - mint_tomorrow) / tkk)
//      For the time from midnight to sunrise, similar equations are used, but
//      the minimum
//  temperature of this day (mint_today) is used instead of mint_tomorrow, and
//  the maximum temperature of the previous day (maxt_yesterday) is used instead
//  of maxt_today. Also, (suns-ti-24) is used for the time variable instead of
//  (suns-ti).
//      These exponential equations for night-time temperature ensure that the
//      curve will
//  be continuous with the daytime equation at sunset, and will pass through the
//  minimum temperature at sunrise.
//
//  Input argument:
//     ti - time of day (hours).
//  Global variables used:
//     DayLength, Daynum, LastDayWeatherData, pi, SitePar, SolarNoon, sunr, suns
//
unsafe fn daytmp(ti: f64) -> f64 {
    const tkk: f64 = 15.; // The temperature increase at which the sensible heat flux is doubled, in comparison with the situation without buoyancy.
    const tcoef: f64 = 4.; // time coefficient for the exponential part of the equation.
    let hmax = SolarNoon + SitePar[8]; // hour of maximum temperature
    let im1 = Daynum - 1; // day of year yesterday
    let mut ip1 = Daynum + 1; // day of year tomorrow
    if ip1 > LastDayWeatherData {
        ip1 = Daynum;
    }
    let amp; // amplitude of temperatures for a period.
    let sst; // the temperature at sunset.
    let st; // computed from time of day, used for daytime temperature.
    let sts; // intermediate variable for computing sst.
    let HourlyTemperature; // computed temperature at time ti.
                           //
    if ti <= sunr {
        //  from midnight to sunrise
        amp = (GetFromClim(CLIMATE_METRIC_TMAX, im1) - GetFromClim(CLIMATE_METRIC_TMIN, Daynum))
            * (1.
                + (GetFromClim(CLIMATE_METRIC_TMAX, im1)
                    - GetFromClim(CLIMATE_METRIC_TMIN, Daynum))
                    / tkk);
        sts = (std::f64::consts::PI * DayLength / (DayLength + 2. * SitePar[8])).sin();
        //  compute temperature at sunset:
        sst = GetFromClim(CLIMATE_METRIC_TMIN, Daynum) - tkk / 2.
            + 0.5 * (tkk * tkk + 4. * amp * tkk * sts).sqrt();
        HourlyTemperature = (GetFromClim(CLIMATE_METRIC_TMIN, Daynum)
            - sst * ((DayLength - 24.) / tcoef).exp()
            + (sst - GetFromClim(CLIMATE_METRIC_TMIN, Daynum)) * ((suns - ti - 24.) / tcoef).exp())
            / (1. - ((DayLength - 24.) / tcoef).exp());
    } else if ti <= hmax {
        //  from sunrise to hmax
        amp = (GetFromClim(CLIMATE_METRIC_TMAX, Daynum) - GetFromClim(CLIMATE_METRIC_TMIN, Daynum))
            * (1.
                + (GetFromClim(CLIMATE_METRIC_TMAX, Daynum)
                    - GetFromClim(CLIMATE_METRIC_TMIN, Daynum))
                    / tkk);
        st = (std::f64::consts::PI * (ti - SolarNoon + DayLength / 2.)
            / (DayLength + 2. * SitePar[8]))
            .sin();
        HourlyTemperature = GetFromClim(CLIMATE_METRIC_TMIN, Daynum) - tkk / 2.
            + 0.5 * (tkk * tkk + 4. * amp * tkk * st).sqrt();
    } else if ti <= suns {
        //  from hmax to sunset
        amp = (GetFromClim(CLIMATE_METRIC_TMAX, Daynum) - GetFromClim(CLIMATE_METRIC_TMIN, ip1))
            * (1.
                + (GetFromClim(CLIMATE_METRIC_TMAX, Daynum)
                    - GetFromClim(CLIMATE_METRIC_TMIN, ip1))
                    / tkk);
        st = (std::f64::consts::PI * (ti - SolarNoon + DayLength / 2.)
            / (DayLength + 2. * SitePar[8]))
            .sin();
        HourlyTemperature = GetFromClim(CLIMATE_METRIC_TMIN, ip1) - tkk / 2.
            + 0.5 * (tkk * tkk + 4. * amp * tkk * st).sqrt();
    } else
    //  from sunset to midnight
    {
        amp = (GetFromClim(CLIMATE_METRIC_TMAX, Daynum) - GetFromClim(CLIMATE_METRIC_TMIN, ip1))
            * (1.
                + (GetFromClim(CLIMATE_METRIC_TMAX, Daynum)
                    - GetFromClim(CLIMATE_METRIC_TMIN, ip1))
                    / tkk);
        sts = (std::f64::consts::PI * DayLength / (DayLength + 2. * SitePar[8])).sin();
        sst = GetFromClim(CLIMATE_METRIC_TMIN, ip1) - tkk / 2.
            + 0.5 * (tkk * tkk + 4. * amp * tkk * sts).sqrt();
        HourlyTemperature = (GetFromClim(CLIMATE_METRIC_TMIN, ip1)
            - sst * ((DayLength - 24.) / tcoef).exp()
            + (sst - GetFromClim(CLIMATE_METRIC_TMIN, ip1)) * ((suns - ti) / tcoef).exp())
            / (1. - ((DayLength - 24.) / tcoef).exp());
    }
    return HourlyTemperature;
    //     Reference:
    //     Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling
    //  diurnal patterns of air temperature, radiation, wind speed and
    //  relative humidity by equations from daily characteristics.
    //  Agricultural Systems 51:377-393.
}

//     Function tdewhour() computes the hourly values of dew point
//  temperature from average dew-point and the daily estimated range.
//  This range is computed as a regression on maximum and minimum temperatures.
//  Input arguments:
//     ti - time of day (hours).
//     tt - air temperature C at this time of day.
//  Global variables used:
//    Daynum, LastDayWeatherData, SitePar, SolarNoon, sunr, suns.
//
unsafe fn tdewhour(ti: f64, tt: f64) -> f64 {
    let im1 = Daynum - 1; // day of year yeaterday
    let mut ip1 = Daynum + 1; // day of year tomorrow
    if ip1 > LastDayWeatherData {
        ip1 = Daynum;
    }
    let tdewhr; // the dew point temperature (c) of this hour.
    let tdmin; // minimum of dew point temperature.
    let mut tdrange; // range of dew point temperature.
    let hmax = SolarNoon + SitePar[8]; // time of maximum air temperature

    if ti <= sunr {
        // from midnight to sunrise
        tdrange = SitePar[12]
            + SitePar[13] * GetFromClim(CLIMATE_METRIC_TMAX, im1)
            + SitePar[14] * GetFromClim(CLIMATE_METRIC_TMIN, Daynum);
        if tdrange < 0. {
            tdrange = 0.;
        }
        tdmin = GetFromClim(CLIMATE_METRIC_TDEW, im1) - tdrange / 2.;
        tdewhr = tdmin
            + tdrange * (tt - GetFromClim(CLIMATE_METRIC_TMIN, Daynum))
                / (GetFromClim(CLIMATE_METRIC_TMAX, im1)
                    - GetFromClim(CLIMATE_METRIC_TMIN, Daynum));
    } else if ti <= hmax {
        // from sunrise to hmax
        tdrange = SitePar[12]
            + SitePar[13] * GetFromClim(CLIMATE_METRIC_TMAX, Daynum)
            + SitePar[14] * GetFromClim(CLIMATE_METRIC_TMIN, Daynum);
        if tdrange < 0. {
            tdrange = 0.;
        }
        tdmin = GetFromClim(CLIMATE_METRIC_TDEW, Daynum) - tdrange / 2.;
        tdewhr = tdmin
            + tdrange * (tt - GetFromClim(CLIMATE_METRIC_TMIN, Daynum))
                / (GetFromClim(CLIMATE_METRIC_TMAX, Daynum)
                    - GetFromClim(CLIMATE_METRIC_TMIN, Daynum));
    } else {
        //  from hmax to midnight
        tdrange = SitePar[12]
            + SitePar[13] * GetFromClim(CLIMATE_METRIC_TMAX, Daynum)
            + SitePar[14] * GetFromClim(CLIMATE_METRIC_TMIN, ip1);
        if tdrange < 0. {
            tdrange = 0.;
        }
        tdmin = GetFromClim(CLIMATE_METRIC_TDEW, ip1) - tdrange / 2.;
        tdewhr = tdmin
            + tdrange * (tt - GetFromClim(CLIMATE_METRIC_TMIN, ip1))
                / (GetFromClim(CLIMATE_METRIC_TMAX, Daynum)
                    - GetFromClim(CLIMATE_METRIC_TMIN, ip1));
    }
    return tdewhr;
}

///     Function dayrh() computes the hourly values of relative humidity, using
//     the hourly air
//  and dew point temperatures. It calls function VaporPressure().
//     If the estimated dew point is higher than the actual air temperature, its
//     value is
//  taken as the air temperature (relative humidity 100%).
//     The relative humidity is calculated as the percentage ratio of the
//     saturated vapor
//  pressure at dew point temperature and the saturated vapor pressure at actual
//  air temperature.
//     Input arguments:
//        tt - air temperature C at this time of day.
//        tdew - dew point temperature C at this time of day.
//
fn dayrh(tt: f64, tdew: f64) -> f64 {
    let td = fmin(tt, tdew); // the dew point temperature (C), is assumed to
                             // be tt if tt < tdew.
    let esvp = VaporPressure(tt); // the saturated vapor pressure in the air (mbar).
    let vpa = VaporPressure(td); // the actual vapor pressure in the air (mbar).
    let mut RelHumHour = 100. * vpa / esvp; // relative humidity at this time of day, %.
    if RelHumHour < 1. {
        RelHumHour = 1.;
    }
    if RelHumHour > 100. {
        RelHumHour = 100.;
    }
    return RelHumHour;
    //     Reference:
    //     Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling
    //  diurnal patterns of air temperature, radiation, wind speed and
    //  relative humidity by equations from daily characteristics.
    //  Agricultural Systems 51:377-393.
}

//     The daywnd function computes the hourly values of wind speed
//  (m / sec), estimated from the measured total daily wind run.
//     Input arguments:
//       t1 = the hour at which day-time wind begins to blow.
//       t2 = the hour at which day-time wind speed is maximum.
//       t3 = the hour at which day-time wind ceases to blow.
//       ti = the hour of the day.
//       wind = the daily total wind run (km per day).
//       wnytf = Factor for estimating night-time wind (from time
//               t3 to time t1 next day).
//
//     The algorithm is described by Ephrath et al. (1996). It is based on the
//	following assumptions:
//     Although the variability of wind speed during any day is very
//  large, the diurnal wind speed curves appear to be characterized by
//  the following repetitive pattern: increase in wind speed from time
//  t1 in the morning to time t2 in the afternoon, decrease from t2 to t3
//  in the evening, and a low constant wind speed at night, from t3 to t1
//  in the next day.
//     The values of t1, t2, and t3 have been determined in the calling routine:
//  t1 is SitePar(1) hours after sunrise, t2 is SitePar(2) hours after solar
//  noon, and t3 is SitePar(3) hours after sunset. These parameters are site-
//  specific. They are 1, 3, and 0, respectively, for the San Joaquin valley of
//  California and for Arizona, and 1, 4, and 2, respectively, for the coastal
//  plain of israel.
//     The wind speed during the night, from t3 to t1 next day (wmin)
//  is assumed to be proportional to the daily total wind run. The
//  ratio wnytf is also site-specific, SitePar(4), ( 0.008 for San Joaquin
//  and Arizona, 0.0025 for the coastal plain of Israel). wmin is the minimum
//  wind speed from t1 to t3.
//      wtday is computed by subtracting the daily integral of wmin, after
//  converting it from m/sec to km/day, from the total daily wind run (wndt).
//      wmax, the maximum wind speed at time t2 (minus wmin), is
//  computed from wtday and converted to m/sec.
//      daywnd from t1 to t2 is now computed as an increasing
//  sinusoidal function from wmin to wmin + wmax, and it is computed from
//  t2 to t3 as a decreasing sinusoidal function from wmin + wmax to  wmin.
//
fn daywnd(ti: f64, wind: f64, t1: f64, t2: f64, t3: f64, wnytf: f64) -> f64 {
    let HourlyWind;
    //   constants related to t1, t2, t3 :
    let sf1 = 4. * (t2 - t1);
    let sf2 = 4. * (t3 - t2);
    let wmin = wind * wnytf; //  the constant minimum wind speed during the night (m/sec).
    let wtday = wind - wmin * 3.6 * 24.; //  integral of wind run from t1 to t3, minus wmin (km).
    let wmax = wtday * 2. * std::f64::consts::PI / 3.6 / (sf1 + sf2); //  the maximum wind speed (m per sec), above wmin.

    if ti >= t1 && ti < t2 {
        HourlyWind = wmin + wmax * (2. * std::f64::consts::PI * (ti - t1) / sf1).sin();
    } else if ti >= t2 && ti < t3 {
        HourlyWind = wmin + wmax * (2. * std::f64::consts::PI * (ti - (2. * t2 - t3)) / sf2).sin();
    } else {
        HourlyWind = wmin;
    }
    return HourlyWind;
    //
    //     Reference:
    //     Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling
    //  diurnal patterns of air temperature, radiation, wind speed and
    //  relative humidity by equations from daily characteristics.
    //  Agricultural Systems 51:377-393.
}

//     Function AverageAirTemperatures() calculates daily average temperatures,
//     daytime
//  average and night time average.
//     Global variables referenced:
//        AirTemp[], Radiation[].
//     Global variables computed:
//        AvrgDailyTemp, NightTimeTemp, DayTimeTemp.
//
unsafe fn AverageAirTemperatures() {
    let mut nn1 = 0; // counter of night hours
    let mut nn2 = 0; // counter of daytime hours
    AvrgDailyTemp = 0.;
    NightTimeTemp = 0.;
    DayTimeTemp = 0.;

    for ihr in 0..24 {
        AvrgDailyTemp += AirTemp[ihr];
        if Radiation[ihr] <= 0. {
            NightTimeTemp += AirTemp[ihr];
            nn1 += 1;
        } else {
            DayTimeTemp += AirTemp[ihr];
            nn2 += 1;
        }
    }
    //
    AvrgDailyTemp = AvrgDailyTemp / 24.;
    NightTimeTemp /= nn1 as f64;
    DayTimeTemp /= nn2 as f64;
}
/// VaporPressure() computes the water vapor pressure in the air (in KPa units) as a function of the air at temperature tt (C). This equation is widely used.
pub fn VaporPressure(tt: f64) -> f64 {
    0.61078 * (17.269 * tt / (tt + 237.3)).exp()
}
///     Function EvapoTranspiration() computes the rate of reference
//     evapotranspiration
//  and related variables. The following subroutines and functions are called
//  for each hour: sunangle, cloudcov(), clcor(), refalbed(), VaporPressure(),
//  clearskyemiss(), del(), gam().
//     Global variables used:
//        declination, Elevation, Latitude, pi, Radiation, RelativeHumidity,
//        SitePar, AirTemp.
//     Global variables computed: CloudCoverRatio, CloudTypeCorr, albedo, rnet,
//     ReferenceETP,
//        es1hour, es2hour, ReferenceTransp.
//     argument: jtout - index indicating if output is required.
unsafe fn EvapoTranspiration() {
    const stefb: f64 = 5.77944E-08; // the Stefan-Boltzman constant, in W m-2 K-4 (= 1.38E-12 * 41880)
    const c12: f64 = 0.125; // c12 ... c15 are constant parameters.
    const c13: f64 = 0.0439;
    const c14: f64 = 0.030;
    const c15: f64 = 0.0576;
    let mut iamhr: usize = 0; // earliest time in day for computing cloud cover
    let mut ipmhr: usize = 0; // latest time in day for computing cloud cover
    let mut cosz: f64 = 0.; // cosine of sun angle from zenith for this hour
    let mut suna: f64 = 0.; // sun angle from horizon, degrees at this hour
    let mut isr: f64; // hourly extraterrestrial radiation in W / m**2

    for ihr in 0..24 as usize {
        let ti = ihr as f64 + 0.5; // middle of the hourly interval
                                   //      The following subroutines and functions are called for each
                                   //  hour: sunangle, cloudcov, clcor, refalbed .
        sunangle(ti, &mut cosz, &mut suna);
        isr = tmpisr * cosz;
        CloudCoverRatio[ihr] = cloudcov(Radiation[ihr], isr, cosz);
        //      clcor is called to compute cloud-type correction.
        //      iamhr and ipmhr are set.
        CloudTypeCorr[ihr] = clcor(ihr, SitePar[7], isr, cosz);
        if (cosz >= 0.1736) && (iamhr == 0) {
            iamhr = ihr;
        }
        if ihr >= 12 && cosz <= 0.1736 && ipmhr == 0 {
            ipmhr = ihr - 1;
        }
        //      refalbed is called to compute the reference albedo for each
        //      hour.
        albedo[ihr] = refalbed(isr, Radiation[ihr], cosz, suna);
    } //   End of 1st hourly loop
      //     Zero some variables that will later be used for summation.
    ReferenceTransp = 0.;
    Rn = 0.; // daily net radiation
    let mut rnet: [f64; 24] = [0.; 24]; // hourly net radiation
    for ihr in 0..24 as usize {
        //      Compute saturated vapor pressure (svp), using function
        //      VaporPressure(). The actual vapor pressure (vp) is computed from
        //      svp and the
        //  relative humidity. Compute vapor pressure deficit (vpd). This
        //  procedure is based on the CIMIS algorithm.
        let svp = VaporPressure(AirTemp[ihr]); // saturated vapor pressure, mb
        let vp = 0.01 * RelativeHumidity[ihr] * svp; // vapor pressure, mb
        let vpd = svp - vp; // vapor pressure deficit, mb.
                            //   Get cloud cover and cloud correction for night hours
        if ihr < iamhr || ihr > ipmhr {
            CloudCoverRatio[ihr] = 0.;
            CloudTypeCorr[ihr] = 0.;
        }
        //     The hourly net radiation is computed using the CIMIS algorithm
        //     (Dong et al., 1988): rlonin, the hourly incoming long wave
        //     radiation, is computed from ea0, cloud cover
        //  (CloudCoverRatio), air temperature (tk),  stefb, and cloud type
        //  correction (CloudTypeCorr).
        //     rnet, the hourly net radiation, W m-2, is computed from the
        //     global radiation, the albedo,
        //  the incoming long wave radiation, and the outgoing longwave
        //  radiation.
        let tk = AirTemp[ihr] + 273.161; // hourly air temperature in Kelvin.
        let ea0 = clearskyemiss(vp, tk); // clear sky emissivity for long wave radiation
                                         //     Compute incoming long wave radiation:
        let rlonin =
            (ea0 * (1. - CloudCoverRatio[ihr]) + CloudCoverRatio[ihr]) * stefb * tk.powi(4)
                - CloudTypeCorr[ihr];
        rnet[ihr] = (1. - albedo[ihr]) * Radiation[ihr] + rlonin - stefb * tk.powi(4);
        Rn += rnet[ihr];
        //     The hourly reference evapotranspiration ReferenceETP is computed
        //     by the
        //  CIMIS algorithm using the modified Penman equation:
        //     The weighting ratio (w) is computed from the functions del() (the
        //     slope of the saturation
        //  vapor pressure versus air temperature) and gam() (the psychometric
        //  constant).
        let w = del(tk, svp) / (del(tk, svp) + gam(Elevation, AirTemp[ihr])); // coefficient of the Penman equation
                                                                              //     The wind function (fu2) is computed using different sets of
                                                                              //     parameters
                                                                              //  for day-time and night-time. The parameter values are as suggested
                                                                              //  by CIMIS.
                                                                              // wind function for computing evapotranspiration
        let fu2 = if Radiation[ihr] <= 0. {
            c12 + c13 * WindSpeed[ihr]
        } else {
            c14 + c15 * WindSpeed[ihr]
        };
        //     hlathr, the latent heat for evaporation of water (W m-2 per mm at
        //     this hour) is
        //  computed as a function of temperature.
        let hlathr = 878.61 - 0.66915 * (AirTemp[ihr] + 273.161);
        //     ReferenceETP, the hourly reference evapotranspiration, is now
        //     computed by the
        //  modified Penman equation.
        ReferenceETP[ihr] = w * rnet[ihr] / hlathr + (1. - w) * vpd * fu2;
        if ReferenceETP[ihr] < 0. {
            ReferenceETP[ihr] = 0.;
        }
        //     ReferenceTransp is the sum of ReferenceETP
        ReferenceTransp += ReferenceETP[ihr];
        //     es1hour and es2hour are computed as the hourly potential
        //     evapotranspiration due to
        //  radiative and aerodynamic factors, respectively.
        //     es1hour and ReferenceTransp are not computed for periods of
        //     negative net radiation.
        es2hour[ihr] = (1. - w) * vpd * fu2;
        if rnet[ihr] > 0. {
            es1hour[ihr] = w * rnet[ihr] / hlathr;
        } else {
            es1hour[ihr] = 0.;
        }
    } //   end of 2nd hourly loop
}

/// Function clearskyemiss() estimates clear sky emissivity for long wave radiation.
/// Input arguments:
/// * `vp` - vapor pressure of the air in KPa
/// * `tk` - air temperature in K.
pub fn clearskyemiss(vp: f64, tk: f64) -> f64 {
    // Convert vp to mbars
    // vapor pressure of the air in mbars.
    let vp1 = vp * 10.;
    // Compute clear sky emissivity by the method of Idso (1981)
    let ea0 = 0.70 + 5.95e-05 * vp1 * (1500. / tk).exp();
    if ea0 > 1. {
        1.
    } else {
        ea0
    }
}
//    Reference:
//      Idso, S.B. 1981. A set of equations for full spectrum and 8- to 14-um
// and 10.5- to 12.5- um thermal radiation from cloudless skies. Water
// Resources Res. 17:295.

/// Computes cloud cover for this hour from radiation data, using the CIMIS algorithm.
/// The return value is cloud cover ratio ( 0 to 1 )
/// Input arguments:
/// * `radihr` - hourly global radiation in W m<sup>-2</sup>.
/// * `isr` - hourly extraterrestrial radiation in W m<sup>-2</sup>.
/// * `cosz` - cosine of sun angle from zenith.
///
/// This algorithm is described by Dong et al. (1988).
/// Cloud cover fraction is estimated as a function of the ratio of actual solar radiation to extraterrestrial radiation.
/// The parameters of this function have been based on California data.
///
/// The equation is for daylight hours, when the sun is not less than 10 degrees above the horizon (coszhr > 0.1736).
fn cloudcov(radihr: f64, isr: f64, cosz: f64) -> f64 {
    const p1: f64 = 1.333; //    p1, p2, p3 are constant parameters.
    const p2: f64 = 1.7778;
    const p3: f64 = 0.294118;
    let mut rasi = 0.; // ratio of radihr to isr.
    let mut clcov = 0.; // computed cloud cover.

    if isr > 0. {
        rasi = radihr / isr;
    }
    if cosz > 0.1736 && rasi <= p1 / p2 {
        clcov = (p1 - p2 * if rasi >= 0.375 { rasi } else { 0.375 }).powf(p3);
        if clcov < 0. {
            clcov = 0.;
        }
    }
    return clcov;
}
//      Reference:
//      Dong, A., Prashar, C.K. and Grattan, S.R. 1988. Estimation of
// daily and hourly net radiation. CIMIS Final Report June 1988, pp.
// 58-79.

///     Function clcor() computes cloud type correction, using the CIMIS
//     algorithm. Input arguments:
//            ck = cloud type correction factor (data for this location).
//            coszhr = cosine of sun angle from zenith.
//            ihr = time of day, hours.
//            isrhr = hourly extraterrestrial radiation in W m-2 .
//     Global variables used: DayLength, Radiation[], SolarNoon, pi
//     Note: This algorithm is described by Dong et al. (1988). ck is the
//  cloud-type correction used in the Monteith equation for estimating
//  net radiation. The value of this correction depends on site and
//  time of year. Regional ck values for California are given by Dong
//  et al. (1988). In the San Joaquin valley of California ck is almost
//  constant from April to October, with an average value of 60. The value
//  of ck is site-dependant, assumed to be constant during the growing season.
//     The daily ck is converted to an hourly value for
//  clear or partly cloudy sky (rasi >= 0.375) and when the sun is at
//  least 10 degrees above the horizon.
//     Evening, night and early morning cloud type correction is temporarily
//  assigned 0. It is later assigned the values of first or
// last non-zero values (in the calling routine).
//
unsafe fn clcor(ihr: usize, ck: f64, isrhr: f64, coszhr: f64) -> f64 {
    let mut rasi = 0.; //  ratio of Radiation to isrhr.
    if isrhr > 0. {
        rasi = Radiation[ihr] / isrhr;
    }
    if coszhr >= 0.1736 && rasi >= 0.375 {
        let angle = std::f64::consts::PI * (ihr as f64 - SolarNoon + 0.5) / DayLength; // hour angle (from solar noon) in radians.
        ck * std::f64::consts::PI / 2. * angle.cos()
    } else {
        0.
    }
}
//      Reference:
//      Dong, A., Prashar, C.K. and Grattan, S.R. 1988. Estimation of
// daily and hourly net radiation. CIMIS Final Report June 1988, pp.
// 58-79.

/// Computes the slope of the saturation vapor pressure (svp, in mb) versus air temperature (tk, in K).
/// This algorithm is the same as used by CIMIS.
fn del(tk: f64, svp: f64) -> f64 {
    let a = 10f64.powf(-0.0304 * tk);
    let b = tk.powi(2);
    let c = 10f64.powf(-1302.88 / tk);
    return (6790.5 - 5.02808 * tk + 4916.8 * a * b + 174209. * c) * svp / b;
}

/// Function gam() computes the psychometric constant at elevation (elev), m above sea level, and air temperature, C (tt).
/// This algorithm is the same as used by CIMIS.
fn gam(elev: f64, tt: f64) -> f64 {
    let bp = 101.3 - 0.01152 * elev + 5.44e-07 * elev.powi(2); //  barometric pressure, KPa, at this elevation.
    return 0.000646 * bp * (1. + 0.000946 * tt);
}
/// Computes the reference crop albedo, using the CIMIS algorithm.
///
/// This algorithm is described by Dong et al. (1988). Albedo is estimated as a function of sun elevation above the horizon (suna) for clear or partly cloudy sky (rasi >= 0.375) and when the sun is at least 10 degrees above the horizon.
///
/// For very cloudy sky, or when solar altitude is below 10 degrees, the following albedo value is assumed: (p4)+ 0.26
///
/// Input arguments:
/// * `isrhr` - hourly extraterrestrial radiation in W m-2 .
/// * `rad` - hourly global radiation in W / m-2 .
/// * `coszhr` - cosine of sun angle from zenith.
/// * `sunahr` - sun angle from horizon, degrees.
fn refalbed(isrhr: f64, rad: f64, coszhr: f64, sunahr: f64) -> f64 {
    const p1: f64 = 0.00158; //  p1 ... p4 are constant parameters.
    const p2: f64 = 0.386;
    const p3: f64 = 0.0188;
    const p4: f64 = 0.26;
    let mut rasi = 0.; //   ratio of rad to isrhr

    if isrhr > 0. {
        rasi = rad / isrhr;
    }
    if coszhr > 0.1736 && rasi >= 0.375 {
        let refalb = p1 * sunahr + p2 * (-p3 * sunahr).exp();
        if refalb > p4 {
            p4
        } else {
            refalb
        }
    } else {
        p4
    }
}
//      Reference:
//      Dong, A., Prashar, C.K. and Grattan, S.R. 1988. Estimation of
// daily and hourly net radiation. CIMIS Final Report June 1988, pp. 58-79.

//     sunangle.cpp : computes sun angle for any time of day.
//     Input argument:
//        ti = time of day, hours.
//      Output arguments:
//        coszhr = cosine of sun angle from zenith for this hour.
//        sunahr = sun angle from horizon, degrees.
//
unsafe fn sunangle(ti: f64, coszhr: &mut f64, sunahr: &mut f64) {
    //      The latitude is converted to radians (xlat).
    let xlat = Latitude * std::f64::consts::PI / 180.; // latitude in radians.
                                                       // amplitude of the sine of the solar height, computed as the product of cosines of latitude and declination angles.
    let cd = xlat.cos() * declination.cos();
    // seasonal offset of the sine of the solar height, computed as the product of sines of latitude and declination angles.
    let sd = xlat.sin() * declination.sin();
    let hrangle = 15. * (ti - SolarNoon) * std::f64::consts::PI / 180.; // hourly angle converted to radians
    *coszhr = sd + cd * hrangle.cos();
    if *coszhr <= 0. {
        *coszhr = 0.;
        *sunahr = 0.;
    } else if *coszhr >= 1. {
        *coszhr = 1.;
        *sunahr = 90.;
    } else {
        *sunahr = (coszhr.acos() * 180. / std::f64::consts::PI - 90.).abs();
    }
}

//     This function is called from DayClim() and is executed on each day with
//     raifall more
//  than 2 mm. It computes the runoff and the retained portion of the rainfall.
//  Note: This function is based on the code of GOSSYM. No changes have been
//  made from the original GOSSYM code (except translation to C++). It has not
//  been validated by actual field measurement.
//     It calculates the portion of rainfall that is lost to runoff, and reduces
//     rainfall to the
//  amount which is actually infiltrated into the soil. It uses the soil
//  conservation service method of estimating runoff.
//     References:
//  - Brady, Nyle C. 1984. The nature and properties of soils, 9th ed. Macmillan
//  Publishing Co.
//  - Schwab, Frevert, Edminster, and Barnes. 1981. Soil and water conservation
//  engineering, 3rd ed. John Wiley & Sons, Inc.
//
//     The following global variables are referenced here:
//  ClayVolumeFraction, Daynum, DayStart, Irrig (structure), NumIrrigations,
//  SandVolumeFraction.
//     The argument used here:  rain = today,s rainfall.
//     The return value:  the amount of water (mm) lost by runoff.
unsafe fn SimulateRunoff(rain: f64) -> f64 {
    static mut bFirst: bool = true; // if this is the first time the function is called.
    static mut iGroup: i32 = 0; // soil group number (by clay and sand in upper soil layer)
    static mut d01: f64 = 0.; // Adjustment of curve number for soil groups A,B,C.
                              //     The following is computed only the first time the function is called.
                              //     Infiltration rate is estimated from the percent sand and percent clay
                              //     in the Ap layer.
                              //  If clay content is greater than 35%, the soil is assumed to have a
                              //  higher runoff potential, if clay content is less than 15% and sand is
                              //  greater than 70%, a lower runoff potential is assumed. Other soils
                              //  (loams) assumed moderate runoff potential. No 'impermeable' (group D)
                              //  soils are assumed.  References: Schwab, Brady.
    if bFirst {
        if SandVolumeFraction[0] > 0.70 && ClayVolumeFraction[0] < 0.15 {
            //     Soil group A = 1, low runoff potential
            iGroup = 1;
            d01 = 1.0;
        } else if ClayVolumeFraction[0] > 0.35 {
            //     Soil group C = 3, high runoff potential
            iGroup = 3;
            d01 = 1.14;
        } else {
            //     Soil group B = 2, moderate runoff potential
            iGroup = 2;
            d01 = 1.09;
        }
        bFirst = false;
    }
    //     Loop to accumulate 5-day antecedent rainfall (mm) which will affect
    //     the soil's ability
    //  to accept new rainfall. This also includes all irrigations.
    let mut i01 = Daynum - 5;
    if i01 < DayStart {
        i01 = DayStart;
    }
    let i02 = Daynum;
    let mut PreviousWetting = 0.; // five day total (before this day) of rain and irrigation, mm
    for Dayn in i01..i02 {
        let mut amtirr = 0.; // mm water applied on this day by irrigation
        for i in 0..NumIrrigations as usize {
            if Dayn == Irrig[i].day {
                amtirr = Irrig[i].amount;
            }
        }
        PreviousWetting += amtirr + GetFromClim(CLIMATE_METRIC_RAIN, Dayn);
    }
    // Adjusting curve number for antecedent rainfall conditions.
    let d02: f64 = if PreviousWetting < 36. {
        //  low moisture, low runoff potential.
        if iGroup == 1 {
            0.71
        } else if iGroup == 2 {
            0.78
        } else if iGroup == 3 {
            0.83
        } else {
            1.
        }
    } else if PreviousWetting > 53. {
        //  wet conditions, high runoff potential.
        if iGroup == 1 {
            1.24
        } else if iGroup == 2 {
            1.15
        } else if iGroup == 3 {
            1.1
        } else {
            1.
        }
    } else {
        //  moderate conditions
        1.
    };
    //  Assuming straight rows, and good cropping practice:
    let crvnum = 78.0 * d01 * d02; // Runoff curve number, adjusted for moisture and soil type.
                                   // maximum potential difference between rainfall and runoff.
    let d03 = 25400. / crvnum - 254.;
    if rain <= 0.2 * d03 {
        0.
    } else {
        (rain - 0.2 * d03).powi(2) / (rain + 0.8 * d03)
    }
}
