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
