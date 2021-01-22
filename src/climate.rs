/// Function dayrad() computes the hourly values of global radiation, in W m-2,
/// using the measured daily total global radiation.
///
/// The algorithm follows the paper of Spitters et al. (1986). It assumes
/// that atmospheric transmission of radiation is lower near the margins of
/// the daylight period, because of an increase in the path length through
/// the atmosphere at lower solar heights. Radiation is therefore assumed to be
/// proportional to sinb * (1 + c11 * sinb), where the value of c11 is set as 0.4 .
///
/// References:
///
/// Spitters, C.J.T., Toussaint, H.A.J.M. and Goudriaan, J. 1986.
/// Separating the diffuse and direct component of global radiation and
/// its implications for modeling canopy photosynthesis. Part I.
/// Components of incoming radiation. Agric. For. Meteorol. 38:217-229.
///
/// Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling
/// diurnal patterns of air temperature, radiation, wind speed and
/// relative humidity by equations from daily characteristics.
/// Agricultural Systems 51:377-393.
#[no_mangle]
extern "C" fn dayrad(_ti: f64, radsum: f64, sinb: f64, c11: f64) -> f64
// Input arguments:
//   ti - time of day (hours) at the middle of this hourly period.
//   radsum - daily radiation integral.
//   sinb - sine of the solar elevation.
//   c11 - constant parameter (0.4).
{
    let hourly_radiation = radsum * sinb * (1f64 + c11 * sinb);
    if hourly_radiation < 0. {
        0.
    } else {
        hourly_radiation
    }
}

/// `VaporPressure` computes the water vapor pressure in the air (in KPa units)
/// as a function of the air at temperature tt (C). This equation is widely used.
#[no_mangle]
extern "C" fn VaporPressure(tt: f64) -> f64 {
    0.61078 * std::f64::consts::E.powf(17.269 * tt / (tt + 237.3))
}

/// Function `dayrh` computes the hourly values of relative humidity, using
/// the hourly air and dew point temperatures. It calls function `VaporPressure`
///
/// If the estimated dew point is higher than the actual air temperature, its
/// value is taken as the air temperature (relative humidity 100%).
///
/// The relative humidity is calculated as the percentage ratio of the
/// saturated vapor pressure at dew point temperature and the saturated vapor
/// pressure at actual air temperature.
///
/// Reference:
///
/// Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling diurnal
/// patterns of air temperature, radiation, wind speed and relative humidity
/// by equations from daily characteristics. Agricultural Systems 51:377-393.
#[no_mangle]
extern "C" fn dayrh(tt: f64, tdew: f64) -> f64
// Input arguments:
//   tt - air temperature C at this time of day.
//   tdew - dew point temperature C at this time of day.
{
    let td = if tt < tdew { tt } else { tdew }; // the dew point temperature (C), is assumed to be tt if tt < tdew.
    let esvp = VaporPressure(tt); // the saturated vapor pressure in the air (mbar).
    let vpa = VaporPressure(td); // the actual vapor pressure in the air (mbar).
    let relative_humidity = 100. * vpa / esvp; // relative humidity at this time of day, %.
    if relative_humidity < 1. {
        1.
    } else if relative_humidity > 100. {
        100.
    } else {
        relative_humidity
    }
}

/// Function refalbed() computes the reference crop albedo, using the CIMIS algorithm.
///
/// This algorithm is described by Dong et al. (1988). Albedo is estimated as a function of sun elevation above the horizon (suna) for clear or partly cloudy sky (rasi >= 0.375) and when the sun is at least 10 degrees above the horizon.
///
/// For very cloudy sky, or when solar altitude is below 10 degrees, the following albedo value is assumed: (p4)+ 0.26
///
///  Reference:
/// Dong, A., Prashar, C.K. and Grattan, S.R. 1988. Estimation of daily and hourly net radiation. CIMIS Final Report June 1988, pp. 58-79.
#[no_mangle]
extern "C" fn refalbed(isrhr: f64, rad: f64, coszhr: f64, sunahr: f64) -> f64
// Input arguments:
//   isrhr = hourly extraterrestrial radiation in W m-2 .
//   rad = hourly global radiation in W / m-2 .
//   coszhr = cosine of sun angle from zenith.
//   sunahr = sun angle from horizon, degrees.
{
    let p1 = 0.00158; //  p1 ... p4 are constant parameters.
    let p2 = 0.386;
    let p3 = 0.0188;
    let p4 = 0.26;
    let rasi = if isrhr > 0. { rad / isrhr } else { 0. }; //   ratio of rad to isrhr
    if coszhr > 0.1736 && rasi >= 0.375 {
        let refalb = p1 * sunahr + p2 * std::f64::consts::E.powf(-p3 * sunahr); // the reference albedo
        if refalb > p4 {
            p4
        } else {
            refalb
        }
    } else {
        p4
    }
}

/// Function del() computes the slope of the saturation vapor pressure (svp, in mb) versus air temperature (tk, in K).
/// This algorithm is the same as used by CIMIS.
#[no_mangle]
extern "C" fn del(tk: f64, svp: f64) -> f64 {
    let a = 10f64.powf(-0.0304 * tk);
    let b = tk.powi(2);
    let c = 10f64.powf(-1302.88 / tk);
    (6790.5 - 5.02808 * tk + 4916.8 * a * b + 174209f64 * c) * svp / b
}

/// Function gam() computes the psychometric constant at elevation (elev), m above sea level, and air temperature, C (tt).
///
/// This algorithm is the same as used by CIMIS.
#[no_mangle]
extern "C" fn gam(elev: f64, tt: f64) -> f64 {
    let bp = 101.3 - 0.01152 * elev + 5.44e-07 * elev.powi(2); //  barometric pressure, KPa, at this elevation.
    0.000646 * bp * (1f64 + 0.000946 * tt)
}

/// Function clearskyemiss() estimates clear sky emissivity for long wave radiation.
///
/// Reference:
///
/// Idso, S.B. 1981. A set of equations for full spectrum and 8- to 14-um and 10.5- to 12.5- um thermal radiation from cloudless skies. Water Resources Res. 17:295.
#[no_mangle]
extern "C" fn clearskyemiss(vp: f64, tk: f64) -> f64
// Input arguments:
//   vp - vapor pressure of the air in KPa
//   tk - air temperature in K.
{
    let vp1 = vp * 10f64; // vapor pressure of the air in mbars.

    let ea0 = 0.70 + 5.95e-05 * vp1 * std::f64::consts::E.powf(1500f64 / tk); // Compute clear sky emissivity by the method of Idso (1981)
    if ea0 > 1f64 {
        1f64
    } else {
        ea0
    }
}

/// Function cloudcov() computes cloud cover for this hour from radiation data, using the CIMIS algorithm. The return value is cloud cover ratio ( 0 to 1 )
///
/// This algorithm is described by Dong et al. (1988). Cloud cover fraction is estimated as a function of the ratio of actual solar radiation to extraterrestrial radiation. The parameters of this function have been based on California data.
///
/// The equation is for daylight hours, when the sun is not less than 10 degrees above the horizon (coszhr > 0.1736).
///
/// Reference:
///
/// Dong, A., Prashar, C.K. and Grattan, S.R. 1988. Estimation of daily and hourly net radiation. CIMIS Final Report June 1988, pp. 58-79.
#[no_mangle]
extern "C" fn cloudcov(radihr: f64, isr: f64, cosz: f64) -> f64
// Input arguments:
//   radihr = hourly global radiation in W m-2 .
//   isr = hourly extraterrestrial radiation in W m-2 . 
//   cosz = cosine of sun angle from zenith.
{
    let p1 = 1.333; //    p1, p2, p3 are constant parameters.
    let p2 = 1.7778;
    let p3 = 0.294118;
    let rasi = if isr > 0f64 { radihr / isr } else { 0f64 }; // ratio of radihr to isr.

    if cosz > 0.1736 && rasi <= p1 / p2 {
        // computed cloud cover.
        let clcov = if rasi >= 0.375 {
            (p1 - p2 * rasi).powf(p3)
        } else {
            (p1 - p2 * 0.375).powf(p3)
        };
        if clcov < 0f64 {
            0f64
        } else {
            clcov
        }
    } else {
        0f64
    }
}

/// The `daywnd` function computes the hourly values of wind speed (`m/sec`), estimated from the measured total daily wind run.
///
/// The algorithm is described by Ephrath et al. (1996). It is based on the following assumptions:
///
/// Although the variability of wind speed during any day is very large, the diurnal wind speed curves appear to be characterized by the following repetitive pattern: increase in wind speed from time `t1` in the morning to time `t2` in the afternoon, decrease from `t2` to `t3` in the evening, and a low constant wind speed at night, from `t3` to `t1` in the next day.
///
/// The values of `t1`, `t2`, and `t3` have been determined in the calling routine: `t1` is `SitePar(1)` hours after sunrise, `t2` is `SitePar(2)` hours after solar noon, and `t3` is `SitePar(3)` hours after sunset. These parameters are site-specific. They are 1, 3, and 0, respectively, for the San Joaquin valley of California and for Arizona, and 1, 4, and 2, respectively, for the coastal plain of israel.
///
/// The wind speed during the night, from `t3` to `t1` next day (`wmin`) is assumed to be proportional to the daily total wind run. The ratio `wnytf` is also site-specific, `SitePar(4)`, ( `0.008` for San Joaquin and Arizona, `0.0025` for the coastal plain of Israel). wmin is the minimum wind speed from `t1` to `t3`.
///
/// `wtday` is computed by subtracting the daily integral of wmin, after converting it from m/sec to `km/day`, from the total daily wind run (wndt).
///
/// `wmax`, the maximum wind speed at time `t2` (minus `wmin`), is computed from wtday and converted to `m/sec`.
///
/// `daywnd` from t1 to t2 is now computed as an increasing sinusoidal function from wmin to `wmin + wmax`, and it is computed from t2 to t3 as a decreasing sinusoidal function from `wmin + wmax` to `wmin`.
///
/// Reference:
///
/// Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling diurnal patterns of air temperature, radiation, wind speed and relative humidity by equations from daily characteristics. Agricultural Systems 51:377-393.
#[no_mangle]
extern "C" fn daywnd(ti: f64, wind: f64, t1: f64, t2: f64, t3: f64, wnytf: f64) -> f64
// Input arguments:
// t1 = the hour at which day-time wind begins to blow.
// t2 = the hour at which day-time wind speed is maximum.
// t3 = the hour at which day-time wind ceases to blow.
// ti = the hour of the day.
// wind = the daily total wind run (km per day).
// wnytf = Factor for estimating night-time wind (from time t3 to time t1 next day).
{
    let pi = 3.14159f64;
    //   constants related to t1, t2, t3 :
    let sf1 = 4f64 * (t2 - t1);
    let sf2 = 4f64 * (t3 - t2);
    let wmin = wind * wnytf; //  the constant minimum wind speed during the night (m/sec).
    let wtday = wind - wmin * 3.6 * 24f64; //  integral of wind run from t1 to t3, minus wmin (km).
    let wmax = wtday * 2f64 * pi / 3.6 / (sf1 + sf2); //  the maximum wind speed (m per sec), above wmin.
    if ti >= t1 && ti < t2 {
        wmin + wmax * (2f64 * pi * (ti - t1) / sf1).sin()
    } else if ti >= t2 && ti < t3 {
        wmin + wmax * (2f64 * pi * (ti - (2f64 * t2 - t3)) / sf2).sin()
    } else {
        wmin
    }
}

#[no_mangle]
extern "C" fn tdewest(maxt: f64, site5: f64, site6: f64) -> f64
//     This function estimates the approximate daily average dewpoint temperature when 
//  it is not available. It is called by ReadClimateData().
//     Global variables referenced: SitePar[5] and SitePar[6]
//     Argument used:  maxt = maximum temperature of this day.
//
{
    if maxt <= 20f64 {
        site5
    } else if maxt >= 40f64 {
        site6
    } else {
        ((40f64 - maxt) * site5 + (maxt - 20f64) * site6) / 20f64
    }
}

#[no_mangle]
extern "C" fn clcor(
    ihr: u8,
    ck: f64,
    isrhr: f64,
    coszhr: f64,
    day_length: f64,
    radiation: f64,
    solar_noon: f64,
) -> f64
//     Function clcor() computes cloud type correction, using the CIMIS algorithm.
//     Input arguments:
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
//      Reference:
//      Dong, A., Prashar, C.K. and Grattan, S.R. 1988. Estimation of
// daily and hourly net radiation. CIMIS Final Report June 1988, pp.
// 58-79.
{
    //  ratio of Radiation to isrhr.
    let pi = 3.14159f64;
    let rasi = if isrhr > 0f64 {
        radiation / isrhr
    } else {
        0f64
    };
    if coszhr >= 0.1736 && rasi >= 0.375 {
        let angle = pi * ((ihr as f64) - solar_noon + 0.5) / day_length; // hour angle (from solar noon) in radians.
        ck * pi / 2f64 * angle.cos()
    } else {
        0f64
    }
}

#[no_mangle]
extern "C" fn AverageAirTemperatures(
    air_temperature: *const f64,
    radiation: *const f64,
    average_air_temperature: &mut f64,
    average_daytime_air_temperature: &mut f64,
    average_nighttime_air_temperature: &mut f64,
) {
    let air_temperature = unsafe { std::slice::from_raw_parts(air_temperature, 24) };
    let radiation = unsafe { std::slice::from_raw_parts(radiation, 24) };
    *average_air_temperature = air_temperature.iter().sum::<f64>() / 24f64;
    *average_daytime_air_temperature = 0f64;
    *average_nighttime_air_temperature = 0f64;
    let mut night_hours = 0u8;
    for zipped in radiation.iter().zip(air_temperature) {
        if *zipped.0 <= 0f64 {
            night_hours += 1;
            *average_nighttime_air_temperature += *zipped.1
        } else {
            *average_daytime_air_temperature += *zipped.1
        }
    }
    *average_nighttime_air_temperature /= night_hours as f64;
    *average_daytime_air_temperature /= (24 - night_hours) as f64;
}
