use super::*;
use chrono::prelude::*;
use std::f64::consts::PI;

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
    0.61078 * (17.269 * tt / (tt + 237.3)).exp()
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
fn refalbed(isrhr: f64, rad: f64, coszhr: f64, sunahr: f64) -> f64
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
        let refalb = p1 * sunahr + p2 * (-p3 * sunahr).exp(); // the reference albedo
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
fn del(tk: f64, svp: f64) -> f64 {
    let a = 10f64.powf(-0.0304 * tk);
    let b = tk.powi(2);
    let c = 10f64.powf(-1302.88 / tk);
    (6790.5 - 5.02808 * tk + 4916.8 * a * b + 174209f64 * c) * svp / b
}

/// Function gam() computes the psychometric constant at elevation (elev), m above sea level, and air temperature, C (tt).
///
/// This algorithm is the same as used by CIMIS.
fn gam(elev: f64, tt: f64) -> f64 {
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

    let ea0 = 0.70 + 5.95e-05 * vp1 * (1500f64 / tk).exp(); // Compute clear sky emissivity by the method of Idso (1981)
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
fn cloudcov(radihr: f64, isr: f64, cosz: f64) -> f64
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
    //   constants related to t1, t2, t3 :
    let sf1 = 4f64 * (t2 - t1);
    let sf2 = 4f64 * (t3 - t2);
    let wmin = wind * wnytf; //  the constant minimum wind speed during the night (m/sec).
    let wtday = wind - wmin * 3.6 * 24f64; //  integral of wind run from t1 to t3, minus wmin (km).
    let wmax = wtday * 2f64 * PI / 3.6 / (sf1 + sf2); //  the maximum wind speed (m per sec), above wmin.
    if ti >= t1 && ti < t2 {
        wmin + wmax * (2f64 * PI * (ti - t1) / sf1).sin()
    } else if ti >= t2 && ti < t3 {
        wmin + wmax * (2f64 * PI * (ti - (2f64 * t2 - t3)) / sf2).sin()
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

/// computes cloud type correction, using the CIMIS algorithm.
fn clcor(
    ihr: u8,
    ck: f64,
    isrhr: f64,
    coszhr: f64,
    day_length: f64,
    radiation: f64,
    solar_noon: f64,
) -> f64
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
    let rasi = if isrhr > 0f64 {
        radiation / isrhr
    } else {
        0f64
    };
    if coszhr >= 0.1736 && rasi >= 0.375 {
        let angle = PI * ((ihr as f64) - solar_noon + 0.5) / day_length; // hour angle (from solar noon) in radians.
        ck * PI / 2f64 * angle.cos()
    } else {
        0f64
    }
}

#[no_mangle]
extern "C" fn ComputeDayLength(
    doy: u32,
    year: i32,
    latitude: f64,
    longitude: f64,
    declination: &mut f64,
    tmpisr: &mut f64,
    solar_noon: &mut f64,
    day_length: &mut f64,
    sunr: &mut f64,
    suns: &mut f64,
)
//     Function ComputeDayLength() computes day length, declination, time of
//  solar noon, and extra-terrestrial radiation for this day. The CIMIS
//  (California Irrigation Management Information System) algorithms are used.
//     Global variables referenced here:
//  iyear, Latitude, Longitude, pi,
//     Global variables set here:
//  DayLength, declination
//
{
    //     Convert day of year to corresponding angle in radians (xday).
    let xday =
        2f64 * PI * (doy - 1) as f64 / NaiveDate::from_ymd(year + 1, 1, 1).pred().ordinal() as f64;
    //     Compute declination angle for this day. The equation used here for computing it
    //  is taken from the CIMIS algorithm.
    *declination = 0.006918 - 0.399912 * xday.cos() + 0.070257 * xday.sin()
        - 0.006758 * (2f64 * xday).cos()
        + 0.000907 * (2f64 * xday).sin()
        - 0.002697 * (3f64 * xday).cos()
        + 0.001480 * (3f64 * xday).sin();
    //     Compute extraterrestrial radiation in W m-2. The 'solar constant' (average
    //  value = 1367 W m-2) is corrected for this day's distance between earth and the sun. The
    //  equation used here is from the CIMIS algorithm, which is based on the work of Iqbal (1983).
    *tmpisr = 1367f64
        * (1.00011
            + 0.034221 * xday.cos()
            + 0.00128 * xday.sin()
            + 0.000719 * (2f64 * xday).cos()
            + 0.000077 * (2f64 * xday).sin());
    //     Time of solar noon (SolarNoon) is computed by the CIMIS algorithm,
    //  using a correction for longitude (f), and the date correction (exday).
    //     It is assumed that the time zone is "geographically correct". For
    //  example, longitude between 22.5 and 37.5 East is in time zone GMT+2,
    //  and longitude between 112.5 and 127.5 West is in time zone GMT-8.
    //     All daily times in the model are computed by this method.
    let exday = (0.000075 + 0.001868 * xday.cos()
        - 0.032077 * xday.sin()
        - 0.014615 * (2f64 * xday).cos()
        - 0.04089 * (2f64 * xday).sin())
        * 12.
        / PI;
    let st = 15f64 * (longitude / 15f64).floor();
    let mut f = (longitude - st) / 15f64;
    if longitude > 0f64 {
        if f > 0.5 {
            f -= 1f64;
        }
    } else if f < -0.5 {
        f += 1f64;
    }

    *solar_noon = 12f64 - f - exday;
    //     Compute day length, by commonly used equations, from latitude and declination of
    //  this day. Times of sunrise and of sunset are computed from solar noon and day length.
    let xlat = latitude * PI / 180f64;
    let mut ht = -xlat.tan() * declination.tan();
    if ht > 1f64 {
        ht = 1f64; //  arctic circle
    } else if ht < -1f64 {
        ht = -1f64;
    }
    *day_length = 2f64 * ht.acos() * 12f64 / PI;
    *sunr = *solar_noon - *day_length / 2f64;
    *suns = *sunr + *day_length;
}

#[test]
fn test_comput_day_length() {
    /*{
        "sunrise": "2020-10-20T00:52:07+00:00",
        "sunset": "2020-10-20T11:46:53+00:00",
        "solar_noon": "2020-10-20T06:19:30+00:00",
        "day_length": 39286,
        "civil_twilight_begin": "2020-10-20T00:24:28+00:00",
        "civil_twilight_end": "2020-10-20T12:14:31+00:00",
        "nautical_twilight_begin": "2020-10-19T23:52:41+00:00",
        "nautical_twilight_end": "2020-10-20T12:46:19+00:00",
        "astronomical_twilight_begin": "2020-10-19T23:21:05+00:00",
        "astronomical_twilight_end": "2020-10-20T13:17:55+00:00",
    }*/
    let mut declination: f64 = 0.;
    let mut tmpisr: f64 = 0.;
    let mut day_length: f64 = 0.;
    let mut sunrise: f64 = 0.;
    let mut sunset: f64 = 0.;
    let mut solar_noon: f64 = 0.;
    let latitude = 40.54778;
    let longitude = 81.29;
    let yo = 294;
    let year = 2020;
    ComputeDayLength(
        yo,
        year,
        latitude,
        longitude,
        &mut declination,
        &mut tmpisr,
        &mut solar_noon,
        &mut day_length,
        &mut sunrise,
        &mut sunset,
    );
    let fixed_offset = FixedOffset::east((longitude.div_euclid(15.).floor() as i32) * 3600);
    let solar_noon = fixed_offset.yo(year, yo).and_hms(
        solar_noon.floor() as u32,
        ((solar_noon - solar_noon.floor()) * 60.).floor() as u32,
        ((solar_noon * 60. - (solar_noon * 60.).floor()) * 60.).floor() as u32,
    );
    let sunrise = fixed_offset.yo(year, yo).and_hms(
        sunrise.floor() as u32,
        ((sunrise - sunrise.floor()) * 60.).floor() as u32,
        ((sunrise * 60. - (sunrise * 60.).floor()) * 60.).floor() as u32,
    );
    assert_eq!(
        (solar_noon.with_timezone(&Utc)
            - "2020-10-20T06:19:30+00:00"
                .parse::<DateTime<Utc>>()
                .expect("Parse Error"))
        .num_seconds()
        .abs(),
        4
    );
    assert_eq!(
        (sunrise.with_timezone(&Utc)
            - "2020-10-20T00:52:07+00:00"
                .parse::<DateTime<Utc>>()
                .expect("Parse Error"))
        .num_seconds()
        .abs(),
        149
    );
    let rio = (-22.9110137, -43.2093727);
    let yo = 32;
    let year = 2021;
    let mut declination: f64 = 0.;
    let mut tmpisr: f64 = 0.;
    let mut day_length: f64 = 0.;
    let mut sunrise: f64 = 0.;
    let mut sunset: f64 = 0.;
    let mut solar_noon: f64 = 0.;
    ComputeDayLength(
        yo,
        year,
        rio.0,
        rio.1,
        &mut declination,
        &mut tmpisr,
        &mut solar_noon,
        &mut day_length,
        &mut sunrise,
        &mut sunset,
    );
    let fixed_offset = FixedOffset::east((rio.1.div_euclid(15.).floor() as i32) * 3600);
    let sunrise = fixed_offset.yo(year, yo).and_hms(
        sunrise.floor() as u32,
        ((sunrise - sunrise.floor()) * 60.).floor() as u32,
        ((sunrise * 60. - (sunrise * 60.).floor()) * 60.).floor() as u32,
    );
    assert_eq!(
        (sunrise.with_timezone(&Utc)
            - "2021-02-01T05:33:24-03:00"
                .parse::<DateTime<Utc>>()
                .expect("Parse Error"))
        .num_seconds()
        .abs(),
        137
    );
    let sunset = fixed_offset.yo(year, yo).and_hms(
        sunset.floor() as u32,
        ((sunset - sunset.floor()) * 60.).floor() as u32,
        ((sunset * 60. - (sunset * 60.).floor()) * 60.).floor() as u32,
    );
    assert_eq!(
        (sunset.with_timezone(&Utc)
            - "2021-02-01T18:39:38-03:00"
                .parse::<DateTime<Utc>>()
                .expect("Parse Error"))
        .num_seconds()
        .abs(),
        198
    );
}

#[no_mangle]
extern "C" fn AverageAirTemperatures(
    hours: &[Hour; 24usize],
    average_air_temperature: &mut f64,
    average_daytime_air_temperature: &mut f64,
    average_nighttime_air_temperature: &mut f64,
) {
    *average_air_temperature = 0f64;
    *average_daytime_air_temperature = 0f64;
    *average_nighttime_air_temperature = 0f64;
    let mut night_hours = 0u8;
    for hour in hours.iter() {
        if hour.radiation <= 0f64 {
            night_hours += 1;
            *average_nighttime_air_temperature += hour.temperature
        } else {
            *average_daytime_air_temperature += hour.temperature
        }
        *average_air_temperature += hour.temperature;
    }
    *average_air_temperature /= 24f64;
    *average_nighttime_air_temperature /= night_hours as f64;
    *average_daytime_air_temperature /= (24 - night_hours) as f64;
}

/// computes and returns the hourly values of air temperature, using the
/// measured daily maximum and minimum.
///
/// The algorithm is described in Ephrath et al. (1996). It is based on the
/// following assumptions:
///
/// 1. The time of minimum daily temperature is at sunrise.
/// 2. The time of maximum daily temperature is SitePar[8] hours after solar
///    noon.
///  
/// Many models assume a sinusoidal curve of the temperature during the day,
/// but actual data deviate from the sinusoidal curve in the following
/// characteristic way: a faster increase right after sunrise, a near plateau
/// maximum during several hours in the middle of the day, and a rather fast
/// decrease by sunset. The physical reason for this is a more efficient mixing
/// of heated air from ground level into the atmospheric boundary layer, driven
/// by strong lapse temperature gradients buoyancy.
///
/// NOTE: **will be used for "power" as in Fortran notation**.
///
/// A first order approximation is
///
///     daytmp = tmin + (tmax-tmin) * st * tkk / (tkk + daytmp - tmin)
///
/// where
///
///     st = sin(pi * (ti - SolarNoon + dayl / 2) / (dayl + 2 * SitePar[8]))
///
/// Since daytmp appears on both sides of the first equation, it can be solved
/// and written explicitly as:
///
///     daytmp = tmin - tkk/2 + 0.5 * sqrt(tkk**2 + 4 * amp * tkk * st)
///
/// where the amplitude of tmin and tmax is calculated as
///
///     amp = (tmax - tmin) * (1 + (tmax - tmin) / tkk)
///
/// This ensures that temperature still passes through tmin and tmax values.
///
/// The value of tkk was determined by calibration as 15.
///
/// This algorithm is used for the period from sunrise to the time of maximum
/// temperature, hmax. A similar algorithm is used for the time from hmax to
/// sunset, but the value of the minimum temperature of the next day
/// (mint_tomorrow) is used instead of mint_today.
///
/// Night air temperature is described by an exponentially declining curve.
///
/// For the time from sunset to mid-night:
///
///     daytmp = (mint_tomorrow - sst * exp((dayl - 24) / tcoef)
///              + (sst - mint_tomorrow) * exp((suns - ti) / tcoef))
///              / (1 - exp((dayl - 24) / tcoef))
///
/// where tcoef is a time coefficient, determined by calibration as 4, sst is
/// the sunset temperature, determined by the daytime equation as:
///
///     sst = mint_tomorrow - tkk / 2 + 0.5 * sqrt(tkk**2 + 4 * amp * tkk * sts)
///
/// where
///
///     sts  = sin(pi * dayl / (dayl + 2 * SitePar[8]))
///     amp = (tmax - mint_tomorrow) * (1 + (tmax - mint_tomorrow) / tkk)
///
/// For the time from midnight to sunrise, similar equations are used, but the
/// minimum temperature of this day (mint_today) is used instead of
/// mint_tomorrow, and the maximum temperature of the previous day
/// (maxt_yesterday) is used instead of maxt_today. Also, (suns-ti-24) is used
/// for the time variable instead of (suns-ti).
///
/// These exponential equations for night-time temperature ensure that the
/// curve will be continuous with the daytime equation at sunset, and will pass
/// through the minimum temperature at sunrise.
///
/// Reference:
///
/// Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling diurnal patterns of air temperature, radiation, wind speed and relative humidity by equations from daily characteristics. Agricultural Systems 51:377-393.
#[no_mangle]
extern "C" fn daytmp(
    sim: &Simulation,
    u: u32,
    ti: f64,
    site8: f64,
    LastDayWeatherData: u32,
    sunr: f64,
    suns: f64,
) -> f64
//  Input argument:
//     ti - time of day (hours).
//  Global variables used:
//     DayLength, LastDayWeatherData, pi, SitePar, SolarNoon, sunr, suns
{
    let states = unsafe {
        std::slice::from_raw_parts_mut(sim.states, (sim.day_finish - sim.day_start + 1) as usize)
    };
    let state = &mut states[u as usize];
    let tkk = 15f64; // The temperature increase at which the sensible heat flux is
                     //  doubled, in comparison with the situation without buoyancy.
    let tcoef = 4f64; // time coefficient for the exponential part of the equation.
    let hmax = state.solar_noon + site8; // hour of maximum temperature
    let im1 = if u > 1 { u - 1 } else { 0 }; // day of year yesterday
    let yesterday = sim.climate[im1 as usize];
    let today = sim.climate[u as usize];
    let ip1 = if u + 1 > LastDayWeatherData { u } else { u + 1 };
    let tomorrow = sim.climate[ip1 as usize];
    //
    let amp: f64; // amplitude of temperatures for a period.
    let sst: f64; // the temperature at sunset.
    let st: f64; // computed from time of day, used for daytime temperature.
    let sts: f64; // intermediate variable for computing sst.
    let HourlyTemperature: f64; // computed temperature at time ti.
                                //
    if ti <= sunr {
        //  from midnight to sunrise
        amp = (yesterday.Tmax - today.Tmin) * (1f64 + (yesterday.Tmax - today.Tmin) / tkk);
        sts = (PI * state.day_length / (state.day_length + 2f64 * site8)).sin();
        //  compute temperature at sunset:
        sst = today.Tmin - tkk / 2f64 + 0.5 * (tkk * tkk + 4f64 * amp * tkk * sts).sqrt();
        HourlyTemperature = (today.Tmin - sst * ((state.day_length - 24f64) / tcoef).exp()
            + (sst - today.Tmin) * ((suns - ti - 24f64) / tcoef).exp())
            / (1f64 - ((state.day_length - 24f64) / tcoef).exp());
    } else if ti <= hmax {
        //  from sunrise to hmax
        amp = (today.Tmax - today.Tmin) * (1f64 + (today.Tmax - today.Tmin) / tkk);
        st = (PI * (ti - state.solar_noon + state.day_length / 2.)
            / (state.day_length + 2f64 * site8))
            .sin();
        HourlyTemperature =
            today.Tmin - tkk / 2f64 + 0.5 * (tkk * tkk + 4f64 * amp * tkk * st).sqrt();
    } else if ti <= suns {
        //  from hmax to sunset
        amp = (today.Tmax - tomorrow.Tmin) * (1f64 + (today.Tmax - tomorrow.Tmin) / tkk);
        st = (PI * (ti - state.solar_noon + state.day_length / 2f64)
            / (state.day_length + 2f64 * site8))
            .sin();
        HourlyTemperature =
            tomorrow.Tmin - tkk / 2f64 + 0.5 * (tkk * tkk + 4f64 * amp * tkk * st).sqrt();
    } else {
        //  from sunset to midnight
        amp = (today.Tmax - tomorrow.Tmin) * (1f64 + (today.Tmax - tomorrow.Tmin) / tkk);
        sts = (PI * state.day_length / (state.day_length + 2f64 * site8)).sin();
        sst = tomorrow.Tmin - tkk / 2f64 + 0.5 * (tkk * tkk + 4f64 * amp * tkk * sts).sqrt();
        HourlyTemperature = (tomorrow.Tmin - sst * ((state.day_length - 24f64) / tcoef).exp()
            + (sst - tomorrow.Tmin) * ((suns - ti) / tcoef).exp())
            / (1. - ((state.day_length - 24f64) / tcoef).exp());
    }
    HourlyTemperature
}

/// computes sun angle for any time of day.
fn sunangle(
    ti: f64,
    latitude: f64,
    declination: f64,
    solar_noon: f64,
    coszhr: &mut f64,
    sunahr: &mut f64,
)
// Input argument:
// ti = time of day, hours.
// Output arguments:
// coszhr = cosine of sun angle from zenith for this hour.
// sunahr = sun angle from horizon, degrees.
{
    // The latitude is converted to radians (xlat).
    let xlat = latitude * PI / 180.;
    // amplitude of the sine of the solar height, computed as
    // the product of cosines of latitude and declination angles.
    let cd = xlat.cos() * declination.cos();
    // seasonal offset of the sine of the solar height, computed
    // as the product of sines of latitude and declination angles.
    let sd = xlat.sin() * declination.sin();
    let hrangle = 15. * (ti - solar_noon) * PI / 180.; // hourly angle converted to radians
    *coszhr = sd + cd * hrangle.cos();
    if *coszhr <= 0f64 {
        *coszhr = 0f64;
        *sunahr = 0f64;
    } else if *coszhr >= 1f64 {
        *coszhr = 1f64;
        *sunahr = 90f64;
    } else {
        *sunahr = (coszhr.acos() * 180f64 / PI - 90f64).abs();
    }
}

/// Function tdewhour() computes the hourly values of dew point temperature from average dew-point and the daily estimated range. This range is computed as a regression on maximum and minimum temperatures.
#[no_mangle]
extern "C" fn tdewhour(
    sim: &Simulation,
    u: u32,
    last_day_has_weather_data: u32,
    time: f64,
    temperature: f64,
    sunrise: f64,
    solar_noon: f64,
    site8: f64,
    site12: f64,
    site13: f64,
    site14: f64,
) -> f64 {
    let im1 = if u > 1 { u - 1 } else { 0 }; // day of year yesterday
    let yesterday = sim.climate[im1 as usize];
    let today = sim.climate[u as usize];
    let mut ip1 = u + 1; // day of year tomorrow
    if ip1 > last_day_has_weather_data {
        ip1 = u;
    }
    let tomorrow = sim.climate[ip1 as usize];
    let tdmin; // minimum of dew point temperature.
    let mut tdrange; // range of dew point temperature.
    let hmax = solar_noon + site8; // time of maximum air temperature
    if time <= sunrise {
        // from midnight to sunrise
        tdrange = site12 + site13 * yesterday.Tmax + site14 * today.Tmin;
        if tdrange < 0f64 {
            tdrange = 0f64;
        }
        tdmin = yesterday.Tdew - tdrange / 2f64;
        tdmin + tdrange * (temperature - today.Tmin) / (yesterday.Tmax - today.Tmin)
    } else if time <= hmax {
        // from sunrise to hmax
        tdrange = site12 + site13 * today.Tmax + site14 * today.Tmin;
        if tdrange < 0f64 {
            tdrange = 0f64;
        }
        tdmin = today.Tdew - tdrange / 2f64;
        tdmin + tdrange * (temperature - today.Tmin) / (today.Tmax - today.Tmin)
    } else {
        //  from hmax to midnight
        tdrange = site12 + site13 * today.Tmax + site14 * tomorrow.Tmin;
        if tdrange < 0f64 {
            tdrange = 0f64;
        }
        tdmin = tomorrow.Tdew - tdrange / 2f64;
        tdmin + tdrange * (temperature - tomorrow.Tmin) / (today.Tmax - tomorrow.Tmin)
    }
}

enum SoilRunoff {
    Low,
    Moderate,
    High,
}

#[no_mangle]
extern "C" fn SimulateRunoff(
    sim: &Simulation,
    u: u32,
    SandVolumeFraction: f64,
    ClayVolumeFraction: f64,
    NumIrrigations: u32,
) -> f64
//     This function is called from DayClim() and is executed on each day with raifall more
//  than 2 mm. It computes the runoff and the retained portion of the rainfall. Note: This
//  function is based on the code of GOSSYM. No changes have been made from the original GOSSYM
//  code (except translation to C++). It has not been validated by actual field measurement.
//     It calculates the portion of rainfall that is lost to runoff, and reduces rainfall to the
//  amount which is actually infiltrated into the soil. It uses the soil conservation service
//  method of estimating runoff.
//     References:
//  - Brady, Nyle C. 1984. The nature and properties of soils, 9th ed. Macmillan Publishing Co.
//  - Schwab, Frevert, Edminster, and Barnes. 1981. Soil and water conservation engineering,
//  3rd ed. John Wiley & Sons, Inc.
//
//     The following global variables are referenced here:
//  ClayVolumeFraction, Irrig (structure), NumIrrigations, SandVolumeFraction.
//     The argument used here:  rain = today,s rainfall.
//     The return value:  the amount of water (mm) lost by runoff.
{
    let iGroup: SoilRunoff;
    let d01: f64; // Adjustment of curve number for soil groups A,B,C.

    // Infiltration rate is estimated from the percent sand and percent clay in the Ap layer.
    // If clay content is greater than 35%, the soil is assumed to have a higher runoff potential,
    // if clay content is less than 15% and sand is greater than 70%, a lower runoff potential is
    // assumed. Other soils (loams) assumed moderate runoff potential. No 'impermeable' (group D)
    // soils are assumed.  References: Schwab, Brady.

    if SandVolumeFraction > 0.70 && ClayVolumeFraction < 0.15 {
        // Soil group A = 1, low runoff potential
        iGroup = SoilRunoff::Low;
        d01 = 1.0;
    } else if ClayVolumeFraction > 0.35 {
        // Soil group C = 3, high runoff potential
        iGroup = SoilRunoff::High;
        d01 = 1.14;
    } else {
        // Soil group B = 2, moderate runoff potential
        iGroup = SoilRunoff::Moderate;
        d01 = 1.09;
    }
    // Loop to accumulate 5-day antecedent rainfall (mm) which will affect the soil's ability
    // to accept new rainfall. This also includes all irrigations.
    let mut i01 = u as i32 - 5;
    if i01 < 0 {
        i01 = 0;
    }
    let mut PreviousWetting = 0f64; // five day total (before this day) of rain and irrigation, mm
    for Dayn in (i01 as usize)..u as usize {
        let mut amtirr = 0f64; // mm water applied on this day by irrigation
        for i in 0..NumIrrigations as usize {
            if (Dayn as i32) == sim.irrigation[i].day {
                amtirr = sim.irrigation[i].amount;
            }
        }
        PreviousWetting += amtirr + sim.climate[Dayn].Rain;
    }
    //
    let d02; // Adjusting curve number for antecedent rainfall conditions.
    if PreviousWetting < 3f64 {
        //  low moisture, low runoff potential.
        d02 = match iGroup {
            SoilRunoff::Low => 0.71,
            SoilRunoff::Moderate => 0.78,
            SoilRunoff::High => 0.83,
        };
    } else if PreviousWetting > 53f64 {
        //  wet conditions, high runoff potential.
        d02 = match iGroup {
            SoilRunoff::Low => 1.24,
            SoilRunoff::Moderate => 1.15,
            SoilRunoff::High => 1.10,
        };
    } else {
        //  moderate conditions
        d02 = 1.00;
    }
    //  Assuming straight rows, and good cropping practice:
    let mut crvnum = 78.0; // Runoff curve number, unadjusted for moisture and soil type.
    crvnum *= d01 * d02; // adjusted curve number
    let d03 = 25400f64 / crvnum - 254f64; // maximum potential difference between rainfall and runoff.
                                          //
    let rain = sim.climate[u as usize].Rain;
    if rain <= 0.2 * d03 {
        0f64
    } else {
        (rain - 0.2 * d03).powi(2) / (rain + 0.8 * d03)
    }
}
/// Function EvapoTranspiration() computes the rate of reference evapotranspiration and related variables.
#[no_mangle]
extern "C" fn EvapoTranspiration(
    state: &mut State,
    latitude: f64,
    elevation: f64,
    declination: f64,
    tmpisr: f64,
    site7: f64,
) {
    const stefb: f64 = 5.77944E-08; // the Stefan-Boltzman constant, in W m-2 K-4 (= 1.38E-12 * 41880)
    const c12: f64 = 0.125; // c12 ... c15 are constant parameters.
    const c13: f64 = 0.0439;
    const c14: f64 = 0.030;
    const c15: f64 = 0.0576;
    let mut iamhr = 0; // earliest time in day for computing cloud cover
    let mut ipmhr = 0; // latest time in day for computing cloud cover
    let mut cosz: f64 = 0f64; // cosine of sun angle from zenith for this hour
    let mut suna: f64 = 0f64; // sun angle from horizon, degrees at this hour
                              //      Start hourly loop
    for (ihr, hour) in state.hours.iter_mut().enumerate() {
        let ti = ihr as f64 + 0.5; // middle of the hourly interval
                                   //      The following subroutines and functions are called for each
                                   //  hour: sunangle, cloudcov, clcor, refalbed .
        sunangle(
            ti,
            latitude,
            declination,
            state.solar_noon,
            &mut cosz,
            &mut suna,
        );
        let isr = tmpisr * cosz; // hourly extraterrestrial radiation in W / m**2
        hour.cloud_cov = cloudcov(hour.radiation, isr, cosz);
        // clcor is called to compute cloud-type correction.
        // iamhr and ipmhr are set.
        hour.cloud_cor = clcor(
            ihr as u8,
            site7,
            isr,
            cosz,
            state.day_length,
            hour.radiation,
            state.solar_noon,
        );
        if (cosz >= 0.1736) && (iamhr == 0) {
            iamhr = ihr;
        }
        if (ihr >= 12) && (cosz <= 0.1736) && (ipmhr == 0) {
            ipmhr = ihr - 1;
        }
        // refalbed is called to compute the reference albedo for each hour.
        hour.albedo = refalbed(isr, hour.radiation, cosz, suna);
    }
    // Zero some variables that will later be used for summation.
    state.evapotranspiration = 0f64;
    state.net_radiation = 0f64; // daily net radiation
    for (ihr, hour) in &mut state.hours.iter_mut().enumerate() {
        //      Compute saturated vapor pressure (svp), using function VaporPressure().
        //      The actual vapor pressure (vp) is computed from svp and the
        //  relative humidity. Compute vapor pressure deficit (vpd). This
        //  procedure is based on the CIMIS algorithm.
        let svp = VaporPressure(hour.temperature); // saturated vapor pressure, mb
        let vp = 0.01 * hour.humidity * svp; // vapor pressure, mb
        let vpd = svp - vp; // vapor pressure deficit, mb.
                            //   Get cloud cover and cloud correction for night hours
        if ihr < iamhr || ihr > ipmhr {
            hour.cloud_cov = 0f64;
            hour.cloud_cor = 0f64;
        }
        //     The hourly net radiation is computed using the CIMIS algorithm (Dong et al., 1988):
        //     rlonin, the hourly incoming long wave radiation, is computed from ea0, cloud cover
        //  (CloudCoverRatio), air temperature (tk),  stefb, and cloud type correction (CloudTypeCorr).
        //     rnet, the hourly net radiation, W m-2, is computed from the global radiation, the albedo,
        //  the incoming long wave radiation, and the outgoing longwave radiation.
        let tk = hour.temperature + 273.161; // hourly air temperature in Kelvin.
        let ea0 = clearskyemiss(vp, tk); // clear sky emissivity for long wave radiation
                                         //     Compute incoming long wave radiation:
        let rlonin =
            (ea0 * (1f64 - hour.cloud_cov) + hour.cloud_cov) * stefb * tk.powi(4) - hour.cloud_cor;
        let rnet = (1f64 - hour.albedo) * hour.radiation + rlonin - stefb * tk.powi(4);
        state.net_radiation += rnet;
        //     The hourly reference evapotranspiration ReferenceETP is computed by the
        //  CIMIS algorithm using the modified Penman equation:
        //     The weighting ratio (w) is computed from the functions del() (the slope of the saturation
        //  vapor pressure versus air temperature) and gam() (the psychometric constant).
        let w = del(tk, svp) / (del(tk, svp) + gam(elevation, hour.temperature)); // coefficient of the Penman equation

        //     The wind function (fu2) is computed using different sets of parameters
        //  for day-time and night-time. The parameter values are as suggested by CIMIS.
        let fu2 = if hour.radiation <= 0f64 {
            c12 + c13 * hour.wind_speed
        } else {
            c14 + c15 * hour.wind_speed
        }; // wind function for computing evapotranspiration

        // hlathr, the latent heat for evaporation of water (W m-2 per mm at this hour) is computed as a function of temperature.
        let hlathr = 878.61 - 0.66915 * (hour.temperature + 273.161);
        // ReferenceETP, the hourly reference evapotranspiration, is now computed by the modified Penman equation.
        hour.ref_et = w * rnet / hlathr + (1f64 - w) * vpd * fu2;
        if hour.ref_et < 0f64 {
            hour.ref_et = 0f64;
        }
        // ReferenceTransp is the sum of ReferenceETP
        state.evapotranspiration += hour.ref_et;
        // es1hour and es2hour are computed as the hourly potential evapotranspiration due to radiative and aerodynamic factors, respectively.
        // es1hour and ReferenceTransp are not computed for periods of negative net radiation.
        hour.et2 = (1f64 - w) * vpd * fu2;
        hour.et1 = if rnet > 0f64 { w * rnet / hlathr } else { 0f64 };
    }
}
