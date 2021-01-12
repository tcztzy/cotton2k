/// This function computes and returns the potential stem growth of cotton plants.
#[no_mangle]
extern "C" fn PotentialStemGrowth(
    stem_dry_weight: f64,
    days_since_emergence: i32,
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
        var12 * (var13 + var14 * days_since_emergence as f64)
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

/// This function simulates the growth in height of the main stem of cotton plants.
/// It is called from PlantGrowth(). It returns the added plant height (cm).
#[no_mangle]
extern "C" fn AddPlantHeight(
    density_factor: f64,
    physiological_days_increment: f64,
    number_of_pre_fruiting_nodes: u32,
    second_fruiting_branch_code: u32,
    age_of_last_pre_fruiting_node: f64,
    age_of_penultimate_pre_fruiting_node: f64,
    average_physiological_age_of_top_three_nodes: f64,
    water_stress_of_stem: f64,
    carbon_stress: f64,
    nitrogen_stress_of_vegetative: f64,
    pixdz: f64,
    days_since_emergence: i32,
    days_since_emergence_to_start_retroactive_plant_map_adjustment: i32,
    number_of_days_for_retroactive_plant_map_adjustment: i32,
    need_adjust: bool,
    adjusted_add_height_rate: f64,
    var19: f64,
    var20: f64,
    var21: f64,
    var22: f64,
    var23: f64,
    var24: f64,
    var25: f64,
    var26: f64,
) -> f64 {
    //     The following constant parameters are used:
    let vhtpar: [f64; 7] = [1.0, 0.27, 0.60, 0.20, 0.10, 0.26, 0.32];
    let mut addz; // daily plant height growth increment, cm.
                  //     Calculate vertical growth of main stem before the square on the second fruiting branch
                  //  has appeared. Added stem height (addz) is a function of the age of the last prefruiting node.
    if second_fruiting_branch_code == 0 {
        addz = vhtpar[0] - vhtpar[1] * age_of_last_pre_fruiting_node;
        if addz > vhtpar[2] {
            addz = vhtpar[2];
        }
        if addz < 0. {
            addz = 0.;
        }
        // It is assumed that the previous prefruiting node is also capable of growth, and its growth (dz2) is added to addz.
        if number_of_pre_fruiting_nodes > 1 {
            // plant height growth increment due to growth of the second node from the top.
            let mut dz2 = var19 - var20 * age_of_penultimate_pre_fruiting_node;
            if dz2 < 0. {
                dz2 = 0.;
            }
            if dz2 > vhtpar[3] {
                dz2 = vhtpar[3];
            }
            addz += dz2;
        }
        // The effect of water stress on stem height at this stage is less than at a later stage (as modified by vhtpar(4)).
        addz *= 1. - vhtpar[4] * (1. - water_stress_of_stem);
    } else {
        // Calculate vertical growth of main stem after the second square has appeared. Added stem height (addz) is a function of the average age of the upper three main stem nodes.
        addz = var21
            + average_physiological_age_of_top_three_nodes
                * (var22 + var23 * average_physiological_age_of_top_three_nodes);
        if average_physiological_age_of_top_three_nodes > (-0.5 * var22 / var23) {
            addz = var24;
        }
        if addz < var24 {
            addz = var24;
        }
        if addz > var25 {
            addz = var25;
        }
        // addz is affected by water, carbohydrate and nitrogen stresses.
        addz *= water_stress_of_stem;
        addz *= 1. - vhtpar[5] * (1. - carbon_stress);
        addz *= 1. - vhtpar[6] * (1. - nitrogen_stress_of_vegetative);
    }
    // The effect of temperature is expressed by DayInc. there are also effects of pix, plant density, and of a variety-specific calibration parameter (VarPar(26)).
    addz *= var26 * pixdz * physiological_days_increment * density_factor;
    // Apply adjustment to addz if plant map data have been read
    let kdadjustend = days_since_emergence_to_start_retroactive_plant_map_adjustment
        + number_of_days_for_retroactive_plant_map_adjustment;
    if days_since_emergence > days_since_emergence_to_start_retroactive_plant_map_adjustment
        && days_since_emergence <= kdadjustend
    {
        if need_adjust {
            addz *= adjusted_add_height_rate;
        }
    }
    return addz;
}
