/// This function calculates the reduction of potential root growth rate in cells with
/// low oxygen content (high water content). It is called from PotentialRootGrowth().
///
/// It has been adapted from GOSSYM, but the critical value of soil moisture potential
/// for root growth reduction (i.e., water logging conditions) has been changed.
#[no_mangle]
extern "C" fn SoilAirOnRootGrowth(psislk: f64, pore_space: f64, vh2oclk: f64) -> f64
//     The following input arguments are used:
//        poreSpace -  value of PoreSpace (v/v) for this layer.
//        psislk -  value of SoilPsi for this cell.
//        vh2oclk - water content (v/v) of this cell
{
    // Constant parameters:
    let p1 = 0f64;
    let p2 = 1f64;
    let p3 = 0.1f64;
    // The following is actually disabled by the choice of the calibration parameters. It
    // may be redefined when more experimental data become available.
    // Reduced root growth when water content is at pore - space saturation
    // (below water table).
    if vh2oclk >= pore_space {
        p3
    }
    // Effect of oxygen deficiency on root growth (the return value).
    else if psislk > p1 {
        p2
    } else {
        1f64
    }
}

/// This function calculates the reduction of potential root growth rate in cells with low nitrate content. It is called from PotentialRootGrowth().
///
/// It has been adapted from GOSSYM. It is assumed that root growth is reduced when nitrate N content falls below a certain level.
///
/// NOTE: This function actually does nothing. It is disabled by the choice of the constant parameters. It may be redefined when more experimental data become available.
#[no_mangle]
#[allow(unused_variables)]
extern "C" fn SoilNitrateOnRootGrowth(vno3clk: f64) -> f64
// The following argument is used:
//   vno3clk - VolNo3NContent value for this cell
{
    1f64
}

/// This function is called from PotentialRootGrowth(), TapRootGrowth() and
/// LateralRootGrowth(). It computes the effects of soil temperature on the rate
/// growth. It is essentially based on the usage of GOSSYM, but relative values
/// are computed here. The computed value returned by this function is between 0 and 1.
///
/// It is assumed that maximum root growth occurs at or above 30 C, and no root growth
/// occurs at or below 13.5 C. A quadratic response to temperature between these limits
/// is assumed.
///
/// The following argument is used:
///
/// t - Soil temperature (C), daily average.
///
/// The parameters used are p1, p2, p3, with the following results:
/// t =      14    16    18    20    22    24    26    28    30
/// trf =  .053  .261  .443  .600  .731  .837  .917  .971  1.00
///
#[no_mangle]
pub extern "C" fn SoilTemOnRootGrowth(t: f64) -> f64 {
    let p1 = -2.12;
    let p2 = 0.2;
    let p3 = -0.0032;

    if t >= 30. {
        1.
    } else {
        let trf = p1 + t * (p2 + p3 * t);
        if trf > 1. {
            1.
        } else if trf < 0. {
            0.
        } else {
            trf
        }
    }
}
