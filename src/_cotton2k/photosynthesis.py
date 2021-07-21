# parameters used to correct photosynthesis for ambient CO2 concentration.
CO2_PARAMETERS = (
    1.0235,
    1.0264,
    1.0285,
    1.0321,
    1.0335,
    1.0353,
    1.0385,
    1.0403,
    1.0431,
    1.0485,
    1.0538,
    1.0595,
    1.0627,
    1.0663,
    1.0716,
    1.0752,
    1.0784,
    1.0823,
    1.0880,
    1.0923,
    1.0968,
    1.1019,
    1.1087,
    1.1172,
    1.1208,
    1.1243,
    1.1311,
    1.1379,
    1.1435,
    1.1490,
    1.1545,
    1.1601,
    1.1656,
    1.1712,
    1.1767,
    1.1823,
    1.1878,
    1.1934,
    1.1990,
    1.2045,
    1.2101,
    1.2156,
    1.2212,
    1.2267,
    1.2323,
)

START_YEAR = 1960
STOP_YEAR = START_YEAR + len(CO2_PARAMETERS) - 1


def ambient_co2_factor(year):
    if year < START_YEAR:
        return 1
    if year <= STOP_YEAR:
        return CO2_PARAMETERS[year - START_YEAR]
    return CO2_PARAMETERS[-1] + 0.004864 * (year - STOP_YEAR)
