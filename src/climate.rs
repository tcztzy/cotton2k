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
