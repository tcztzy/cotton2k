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
