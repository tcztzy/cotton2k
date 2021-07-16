def temperature_on_leaf_growth_rate(t):
    ra = (
        -1.14277 + t * (0.0910026 - t * 0.00152344)
        if t > 24
        else -0.317136 + t * (0.0300712 - t * 0.000416356)
    )
    return 0 if ra < 0 else ra / 0.2162044


def leaf_resistance_for_transpiration(age: float) -> float:
    """This function computes and returns the resistance of leaves of cotton plants to
    transpiration.

    It is assumed to be a function of leaf age.

    :param age: leaf age in physiological days.
    """
    # The following constant parameters are used:
    afac: float = 160.0  # factor used for computing leaf resistance.
    agehi: float = 94.0  # higher limit for leaf age.
    agelo: float = 48.0  # lower limit for leaf age.
    rlmin: float = 0.5  # minimum leaf resistance.

    if age <= agelo:
        return rlmin
    if age >= agehi:
        return rlmin + (agehi - agelo) * (agehi - agelo) / afac
    ax: float = 2.0 * agehi - agelo  # intermediate variable
    return rlmin + (age - agelo) * (ax - age) / afac
