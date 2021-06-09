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
