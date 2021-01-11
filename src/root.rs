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
