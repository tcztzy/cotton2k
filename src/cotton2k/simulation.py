# pylint: disable=no-name-in-module
from math import exp
from typing import Any

from _cotton2k.simulation import Simulation as CySimulation
from _cotton2k.simulation import SimulationEnd
from _cotton2k.simulation import State as CyState


class State:
    _: CyState

    def __init__(self, state: CyState) -> None:
        self._ = state

    def __getattr__(self, name: str) -> Any:
        try:
            return getattr(self._, name)
        except AttributeError:
            return getattr(self, name)

    def __setattr__(self, name: str, value: Any) -> None:
        if name == "_":
            object.__setattr__(self, name, value)
        else:
            try:
                setattr(self._, name, value)
            except AttributeError:
                setattr(self, name, value)

    def __getitem__(self, name: str) -> Any:
        return getattr(self, name)

    @property
    def plant_weight(self):
        return (
            self.root_weight
            + self.stem_weight
            + self.green_bolls_weight
            + self.green_bolls_burr_weight
            + self.leaf_weight
            + self.petiole_weight
            + self.square_weight
            + self.open_bolls_weight
            + self.open_bolls_burr_weight
            + self.reserve_carbohydrate
        )

    # pylint: disable=too-many-arguments
    def column_shading(
        self,
        row_space,
        plant_row_column,
        column_width,
        max_leaf_area_index,
        relative_radiation_received_by_a_soil_column,
    ):
        zint = 1.0756 * self.plant_height / row_space
        for k in range(20):
            sw = (k + 1) * column_width
            if k <= plant_row_column:
                k0 = plant_row_column - k
                sw1 = sw - column_width / 2
            else:
                sw1 = sw - column_width / 2 - (plant_row_column + 1) * column_width
                k0 = k
            shade = 0
            if sw1 < self.plant_height:
                shade = 1 - (sw1 / self.plant_height) ** 2
                if (
                    self.light_interception < zint
                    and self.leaf_area_index < max_leaf_area_index
                ):
                    shade *= self.light_interception / zint
            relative_radiation_received_by_a_soil_column[k0] = max(0.05, 1 - shade)


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
        dayfd += min(max((hour.temperature - p1) / p2, 0), p3)
    return dayfd / 24.0


def compute_light_interception(
    leaf_area_index: float,
    max_leaf_area_index: float,
    plant_height: float,
    row_space: float,
    /,
    *,
    version: int = 0x0400,
):
    if version < 0x0500:
        zint = 1.0756 * plant_height / row_space
        lfint = (
            0.80 * leaf_area_index
            if leaf_area_index <= 0.5
            else 1 - exp(0.07 - 1.16 * leaf_area_index)
        )
        if lfint > zint:
            light_interception = (zint + lfint) / 2
        elif leaf_area_index < max_leaf_area_index:
            light_interception = lfint
        else:
            light_interception = zint
        return light_interception if light_interception < 1 else 1
    param = max(1.16, -0.1 * plant_height + 8)
    return 1 - exp(-param * leaf_area_index)


class Simulation(CySimulation):
    @property
    def states(self):
        return [
            self.state(i) for i in range((self.stop_date - self.start_date).days + 1)
        ]

    def state(self, i):
        return State(self._state(i))

    def read_input(self, *args, **kwargs):  # pylint: disable=unused-argument
        # pylint: disable=attribute-defined-outside-init
        self._current_state = self._state(0)
        self._init_state()
        super().read_input(*args, **kwargs)
        self._initialize_root_data()

    def run(self):
        try:
            self._simulate()
        except RuntimeError:
            pass

    def _simulate(self):
        days = (self.stop_date - self.start_date).days
        for i in range(days):
            # pylint: disable=attribute-defined-outside-init
            self._current_state = self._state(i)
            self._simulate_this_day(i)
            self._copy_state(i)
        self._simulate_this_day(days)

    # pylint: disable=attribute-defined-outside-init
    def _simulate_this_day(self, u):
        state = State(self._current_state)
        if state.date >= self.emerge_date:
            state.kday = (state.date - self.emerge_date).days + 1
            # pylint: disable=access-member-before-definition
            if state.leaf_area_index > self.max_leaf_area_index:
                self.max_leaf_area_index = state.leaf_area_index
            state.light_interception = compute_light_interception(
                state.leaf_area_index,
                self.max_leaf_area_index,
                state.plant_height,
                self.row_space,
                version=self.version,
            )
            state.column_shading(
                self.row_space,
                self.plant_row_column,
                self._column_width,
                self.max_leaf_area_index,
                self.relative_radiation_received_by_a_soil_column,
            )
        else:
            state.kday = 0
            state.light_interception = 0
            self.relative_radiation_received_by_a_soil_column[:] = 1
        # The following functions are executed each day (also before emergence).
        self._daily_climate(u)  # computes climate variables for today.
        self._soil_temperature(
            u
        )  # executes all modules of soil and canopy temperature.
        self._soil_procedures(u)  # executes all other soil processes.
        self._soil_nitrogen(u)  # computes nitrogen transformations in the soil.
        # The following is executed each day after plant emergence:
        if (
            state.date >= self.emerge_date
            # pylint: disable=access-member-before-definition
            and self.emerge_switch > 0
        ):
            # If this day is after emergence, assign to emerge_switch the value of 2.
            self.emerge_switch = 2
            state.day_inc = physiological_age(
                state.hours
            )  # physiological days increment for this day. computes physiological age
            self._defoliate(u)  # effects of defoliants applied.
            self._stress(u)  # computes water stress factors.
            self._get_net_photosynthesis(u)  # computes net photosynthesis.
            self._growth(u)  # executes all modules of plant growth.
            self._phenology(u)  # executes all modules of plant phenology.
            state.plant_nitrogen(
                self.emerge_date,
                self.state(
                    max(u - 32, (self.emerge_date - self.start_date).days)
                ).stem_weight,
            )  # computes plant nitrogen allocation.
        # Check if the date to stop simulation has been reached, or if this is the last
        # day with available weather data. Simulation will also stop when no leaves
        # remain on the plant.
        if state.kday > 10 and state.leaf_area_index < 0.0002:
            raise SimulationEnd

    @property
    def _column_width(self):
        return self.row_space / 20
