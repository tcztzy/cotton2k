def temperature_on_leaf_growth_rate(t):
    ra = (
        -1.14277 + t * (0.0910026 - t * 0.00152344)
        if t > 24
        else -0.317136 + t * (0.0300712 - t * 0.000416356)
    )
    return 0 if ra < 0 else ra / 0.2162044
