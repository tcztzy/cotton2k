/// This is the temperature function for leaf growth rate. It is called by function `PotentialLeafGrowth`. It is based on the original code of GOSSYM, and the parameters are the same. The argument t is air temperature (°C).
///
/// `ra` is divided by the maximum value. Thus, the function returns values between `0.` and `1.`
///
/// This will result :
///
///     maximum value of TemperatureOnLeafGrowthRate = 1.00  at  t = 29.86747
///         for t = 24    TemperatureOnLeafGrowthRate = 0.766
///         for t = 27    TemperatureOnLeafGrowthRate = 0.953
///         for t = 30    TemperatureOnLeafGrowthRate = 0.999
///         for t = 36    TemperatureOnLeafGrowthRate = 0.737
///         for t = 42    TemperatureOnLeafGrowthRate = 0.0
///      and for t <= 24 :
///         for t = 12    TemperatureOnLeafGrowthRate = 0.0
///         for t = 16    TemperatureOnLeafGrowthRate = 0.269
///         for t = 20    TemperatureOnLeafGrowthRate = 0.549
///         for t = 24    TemperatureOnLeafGrowthRate = 0.768
#[no_mangle]
extern "C" fn TemperatureOnLeafGrowthRate(t: f64) -> f64 {
    // Constant parameters used:
    let par: [f64; 8] = [
        24.,
        -1.14277,
        0.0910026,
        0.00152344,
        -0.317136,
        0.0300712,
        0.000416356,
        0.2162044,
    ];
    // intermediate value for computing temperatureOnLeafGrowthRate.
    let ra = if t > par[0] {
        par[1] + t * (par[2] - t * par[3])
    } else {
        par[4] + t * (par[5] - t * par[6])
    };
    if ra < 0. {
        0.
    } else {
        ra / par[7]
    }
}

/// This function computes and returns the resistance of leaves of cotton
/// plants to transpiration. It is assumed to be a function of leaf age.
/// It is called from LeafWaterPotential.
/// 
/// The input argument (age) is leaf age in physiological days.
#[no_mangle]
extern "C" fn LeafResistance(age: f64) -> f64 {
    // The following constant parameters are used:
    let afac: f64 = 160.; // factor used for computing leaf resistance.
    let agehi: f64 = 94.; // higher limit for leaf age.
    let agelo: f64 = 48.; // lower limit for leaf age.
    let rlmin: f64 = 0.5; // minimum leaf resistance.

    if age <= agelo {
        rlmin
    } else if age >= agehi {
        rlmin + (agehi - agelo) * (agehi - agelo) / afac
    } else {
        let ax: f64 = 2. * agehi - agelo; // intermediate variable
        rlmin + (age - agelo) * (ax - age) / afac
    }
}
