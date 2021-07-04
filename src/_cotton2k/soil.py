def compute_soil_surface_albedo(
    water_content: float,
    field_capacity: float,
    residual_water_content: float,
    upper_albedo: float,
    lower_albedo: float,
) -> float:
    """Computes albedo of the soil surface
    Less soil water content, higher albedo
    :param water_content:
    :type water_content: float
    :param field_capacity:
    :type field_capacity: float
    :param residual_water_content:
    :type residual_water_content: float
    :param upper_albedo:
    :type upper_albedo: float
    :param lower_albedo:
    :type lower_albedo: float
    """
    if water_content <= residual_water_content:
        soil_surface_albedo = upper_albedo
    elif water_content >= field_capacity:
        soil_surface_albedo = lower_albedo
    else:
        soil_surface_albedo = lower_albedo + (upper_albedo - lower_albedo) * (
            field_capacity - water_content
        ) / (field_capacity - residual_water_content)
    return soil_surface_albedo


def compute_incoming_short_wave_radiation(
    radiation: float, intercepted_short_wave_radiation: float, albedo: float
) -> tuple[float, float, float]:
    """SHORT WAVE RADIATION ENERGY BALANCE
    :return: short wave (global) radiation (ly / sec), global radiation absorbed
             by soil surface, global radiation reflected up to the vegetation
    :rtype: tuple[float, float, float]
    """
    # Division by 41880 (= 698 * 60) converts from Joules per sq m to
    # langley (= calories per sq cm) Or: from Watt per sq m to langley per sec.
    rzero = radiation / 41880  # short wave (global) radiation (ly / sec).
    rss0 = rzero * (
        1 - intercepted_short_wave_radiation
    )  # global radiation after passing through canopy
    return rzero, rss0 * (1 - albedo), rss0 * albedo


def root_psi(soil_psi: float) -> float:
    """This function returns the effect of soil moisture in cell_{l,k} on cotton root
    potential growth rate. It is called from PotentialRootGrowth() and uses the matric
    potential of this cell.

    It is assumed that almost no root growth occurs when the soil is dryer than -p1
    (-20 bars), and root growth rate is maximum at a matric potential of -4 bars
    (p2 - p1) or wetter.

    Effect of soil moisture on root growth (the return value).

    root_psi is computed here as an empirical third degree function, with values
    between 0.02 and 1.
    """
    return min(max(((20.0 + soil_psi) / 16.0) ** 3, 0.02), 1)
