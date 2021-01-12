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

/// This function returns the effect of soil  moisture in cell l,k on cotton root potential
/// growth rate. It is called from PotentialRootGrowth() and uses the matric potential of this cell.
#[no_mangle]
extern "C" fn SoilWaterOnRootGrowth(psislk: f64) -> f64
// The following argument is used:
//   psislk - soil water potential (bars) of this cell.
{
    // It is assumed that almost no root growth occurs when the soil is dryer than -p1 (-20 bars), and root growth rate is maximum at a matric potential of -4 bars (p2 - p1) or wetter.
    // effect of soil moisture on root growth (the return value). smf is computed here as an empirical third degree function, with values between 0.02 and 1.
    let smf = ((20. + psislk) / 16.).powi(3);
    if smf < 0.02 {
        0.02
    } else if smf > 1f64 {
        1f64
    } else {
        smf
    }
}
