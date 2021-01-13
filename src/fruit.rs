/// This function computes the effect of air temperature (t) on growth rate of bolls in cotton plants. It is called from `PotentialFruitGrowth`.
///
/// Some values computed by this function:
///
/// | t (°C) | tfr             |
/// |:-------|:----------------|
/// | 12     | 0.              |
/// | 15     | 0.336           |
/// | 20     | 0.751           |
/// | 25     | 0.978           |
/// | 26     | 1.              |
/// | 28.5   | 1.024 (maximum) |
/// | 30     | 1.016           |
/// | 35     | 0.866           |
/// | 40     | 0.527           |
/// | 45     | 0.              |
#[no_mangle]
extern "C" fn TemperatureOnFruitGrowthRate(t: f64) -> f64 {
    let p1 = -2.041;
    let p2 = 0.215;
    let p3 = 0.00377;

    let tfr = p1 + t * (p2 - p3 * t);
    if tfr < 0. {
        0.
    } else {
        tfr
    }
}
