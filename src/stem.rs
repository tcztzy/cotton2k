/// This function computes and returns the potential stem growth of cotton plants.
#[no_mangle]
extern "C" fn PotentialStemGrowth(
    stem_dry_weight: f64,
    days_since_emerge: i32,
    third_fruiting_branch_code: u32,
    density_factor: f64,
    var12: f64,
    var13: f64,
    var14: f64,
    var15: f64,
    var16: f64,
    var17: f64,
    var18: f64,
) -> f64 {
    // There are two periods for computation of potential stem growth:
    // (1) Before the appearance of a square on the third fruiting branch.
    // Potential stem growth is a function of plant age (days from emergence).
    if third_fruiting_branch_code == 0 {
        var12 * (var13 + var14 * days_since_emerge as f64)
    }
    // (2) After the appearance of a square on the third fruiting branch.
    // It is assumed that all stem tissue that is more than 32 days old is not active.
    // Potential stem growth is a function of active stem tissue weight, and plant density.
    else {
        // effect of plant density on stem growth rate.
        let mut _density_factor = 1. - var15 * (1. - density_factor);
        if _density_factor < 0.2 {
            _density_factor = 0.2;
        }
        _density_factor * var16 * (var17 + var18 * stem_dry_weight)
    }
}
