class DaysToFirstSquare:  # pylint: disable=too-few-public-methods
    """This class is tricky for speed"""

    accumulated_temperature = 0
    days = 0
    accumulated_stress = 0

    def __call__(
        self, temperature, water_stress, nitrogen_stress, calibration_parameter
    ):
        # average temperature from day of emergence.
        self.days += 1
        self.accumulated_temperature += temperature
        average_temperature = min(self.accumulated_temperature / self.days, 34)
        # cumulative effect of water and N stresses on date of first square.
        self.accumulated_stress += (
            0.08 * (1 - water_stress) * 0.3 * (1 - nitrogen_stress)
        )

        return (
            132.2 + average_temperature * (-7 + average_temperature * 0.125)
        ) * calibration_parameter - self.accumulated_stress


days_to_first_square = DaysToFirstSquare()


def physiological_age(hours) -> float:
    """computes physiological age

    This function returns the daily 'physiological age' increment, based on hourly
    temperatures. It is called each day by `SimulateThisDay`.
    """
    # The threshold value is assumed to be 12 C (p1). One physiological day is
    # equivalent to a day with an average temperature of 26 C, and therefore the heat
    # units are divided by 14 (p2).

    # A linear relationship is assumed between temperature and heat unit accumulation
    # in the range of 12 C (p1) to 33 C (p2*p3+p1). the effect of temperatures higher
    # than 33 C is assumed to be equivalent to that of 33 C.

    # The following constant Parameters are used in this function:
    p1 = 12.0  # threshold temperature, C
    p2 = 14.0  # temperature, C, above p1, for one physiological day.
    p3 = 1.5  # maximum value of a physiological day.

    dayfd = 0.0  # the daily contribution to physiological age (return value).
    for hour in hours:
        # add the hourly contribution to physiological age.
        dayfd += min(max((hour["temperature"] - p1) / p2, 0), p3)
    return dayfd / 24.0
