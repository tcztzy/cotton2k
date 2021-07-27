# pylint: disable=no-name-in-module
import datetime
from enum import IntEnum
from functools import cached_property
from typing import Any

from _cotton2k.simulation import Simulation as CySimulation
from _cotton2k.simulation import SimulationEnd
from _cotton2k.simulation import State as CyState
from cotton2k.photo import Photosynthesis


class Stage(IntEnum):
    NotYetFormed = 0
    Square = 1
    GreenBoll = 2
    MatureBoll = 3
    AbscisedAsBoll = 4
    AbscisedAsSquare = 5
    AbscisedAsFlower = 6
    YoungGreenBoll = 7


class State(Photosynthesis):  # pylint: disable=too-many-instance-attributes
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
            + self.stem_weight  # pylint: disable=no-member
            + self.green_bolls_weight
            + self.green_bolls_burr_weight
            + self.leaf_weight
            + self.petiole_weight
            + self.square_weight
            + self.open_bolls_weight
            + self.open_bolls_burr_weight
            + self.reserve_carbohydrate
        )

    @property
    def agetop(self):
        l = len(self.vegetative_branches[0].fruiting_branches) - 1
        # average physiological age of top three nodes.
        if l < 0:
            return 0
        if l == 0:
            return self.vegetative_branches[0].fruiting_branches[0].nodes[0].age
        if l == 1:
            return (
                self.vegetative_branches[0].fruiting_branches[0].nodes[0].age * 2
                + self.vegetative_branches[0].fruiting_branches[1].nodes[0].age
            ) / 3
        return (
            self.vegetative_branches[0].fruiting_branches[l].nodes[0].age
            + self.vegetative_branches[0].fruiting_branches[l - 1].nodes[0].age
            + self.vegetative_branches[0].fruiting_branches[l - 2].nodes[0].age
        ) / 3


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


class Simulation(CySimulation):
    @cached_property
    def states(self):
        return [
            self.state(i) for i in range((self.stop_date - self.start_date).days + 1)
        ]

    def state(self, i):
        if isinstance(i, datetime.date):
            i = (i - self.start_date).days
        return State(self._state(i))

    def _copy_state(self, i):
        super()._copy_state(i)
        pre = self._state(i)
        post = self._state(i + 1)
        for attr in (
            "carbon_allocated_for_root_growth",
            "delay_of_emergence",
            "extra_carbon",
            "fiber_length",
            "fiber_strength",
            "hypocotyl_length",
            "lint_yield",
            "net_photosynthesis",
            "pavail",
            "seed_layer_number",
            "seed_moisture",
        ):
            setattr(post, attr, getattr(pre, attr))

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

    # pylint: disable=attribute-defined-outside-init,no-member
    def _simulate_this_day(self, u):
        state = State(self._current_state)
        if state.date >= self.emerge_date:
            state.kday = (state.date - self.emerge_date).days + 1
            # pylint: disable=access-member-before-definition
            if state.leaf_area_index > self.max_leaf_area_index:
                self.max_leaf_area_index = state.leaf_area_index
            state.light_interception = state.compute_light_interception(
                self.max_leaf_area_index,
                self.row_space,
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
            old_stem_days = 32
            if state.kday > old_stem_days:
                old_stem_weight = self.state(u - 32).stem_weight
                new_stem_weight = state.stem_weight - old_stem_weight
                growing_stem_weight = new_stem_weight
            else:
                old_stem_weight = 0
                new_stem_weight = (
                    state.stem_weight - self.state(self.emerge_date).stem_weight
                )
                growing_stem_weight = state.stem_weight
            state.get_net_photosynthesis(
                self.climate[u]["Rad"],
                self.per_plant_area,
                self.ptsred,
                old_stem_weight,
            )  # computes net photosynthesis.
            self._growth(u, new_stem_weight)  # executes all modules of plant growth.
            self._phenology(u)  # executes all modules of plant phenology.
            state.plant_nitrogen(
                self.emerge_date, growing_stem_weight
            )  # computes plant nitrogen allocation.
        # Check if the date to stop simulation has been reached, or if this is the last
        # day with available weather data. Simulation will also stop when no leaves
        # remain on the plant.
        if state.kday > 10 and state.leaf_area_index < 0.0002:
            raise SimulationEnd

    def _growth(self, u, new_stem_weight):
        state = State(self._current_state)
        # Call _potential_leaf_growth() to compute potential growth rate of leaves.
        self._potential_leaf_growth(u)
        # If it is after first square, call _potential_fruit_growth() to compute
        # potential growth rate of squares and bolls.
        if (
            len(state.vegetative_branches[0].fruiting_branches) > 0
            and len(state.vegetative_branches[0].fruiting_branches[0].nodes) > 0
            and state.vegetative_branches[0].fruiting_branches[0].nodes[0].stage
            != Stage.NotYetFormed
        ):
            self._potential_fruit_growth(u)
        # Call PotentialStemGrowth() to compute PotGroStem, potential growth rate of
        # stems.
        # The effect of temperature is introduced, by multiplying potential growth rate
        # by day_inc.
        # Stem growth is also affected by water stress(water_stress_stem).
        # stem_potential_growth is limited by (maxstmgr * per_plant_area) g per plant
        # per day.
        maxstmgr = 0.067  # maximum posible potential stem growth, g dm - 2 day - 1.
        state.stem_potential_growth = min(
            maxstmgr * self.per_plant_area,
            state.potential_stem_growth(
                new_stem_weight, self.density_factor, *self.cultivar_parameters[12:19]
            )
            * state.day_inc
            * state.water_stress_stem,
        )
        # Call PotentialRootGrowth() to compute potential growth rate on roots.
        # total potential growth rate of roots in g per slab.this is computed in
        # potential_root_growth() and used in actual_root_growth().
        sumpdr = state.potential_root_growth(
            3, state.soil.number_of_layers_with_root, self.per_plant_area
        )
        # Total potential growth rate of roots is converted from g per slab(sumpdr) to
        # g per plant (state.root_potential_growth).
        # Limit state.root_potential_growth to(maxrtgr * per_plant_area) g per plant
        # per day.
        maxrtgr = 0.045  # maximum possible potential root growth, g dm - 2 day - 1.
        state.root_potential_growth = min(
            maxrtgr * self.per_plant_area,
            sumpdr * 100 * self.per_plant_area / self.row_space,
        )
        # Call dry_matter_balance() to compute carbon balance, allocation of carbon to
        # plant parts, and carbon stress.
        vratio = state.dry_matter_balance(self.per_plant_area)
        # If it is after first square, call actual_fruit_growth() to compute actual
        # growth rate of squares and bolls.
        if (
            len(state.vegetative_branches[0].fruiting_branches) > 0
            and len(state.vegetative_branches[0].fruiting_branches[0].nodes) > 0
            and state.vegetative_branches[0].fruiting_branches[0].nodes[0].stage
            != Stage.NotYetFormed
        ):
            state.actual_fruit_growth()
        # Initialize state.leaf_weight.It is assumed that cotyledons fall off at time
        # of first square. Also initialize state.leaf_area and state.petiole_weight.
        if self.first_square_date is not None:
            state.leaf_weight = 0
            state.leaf_area = 0
        else:
            cotylwt = 0.20  # weight of cotyledons dry matter.
            state.leaf_weight = cotylwt
            state.leaf_area = 0.6 * cotylwt
        state.petiole_weight = 0
        # Call actual_leaf_growth to compute actual growth rate of leaves and compute
        # leaf area index.
        state.actual_leaf_growth(vratio)
        state.leaf_area_index = state.leaf_area / self.per_plant_area
        # Add actual_stem_growth to state.stem_weight.
        state.stem_weight += state.actual_stem_growth
        # Plant density affects growth in height of tall plants.
        htdenf = (
            55  # minimum plant height for plant density affecting growth in height.
        )
        z1 = min(
            max((state.plant_height - htdenf) / htdenf, 0),
            1,
        )  # intermediate variable to compute denf2.
        denf2 = 1 + z1 * (
            self.density_factor - 1
        )  # effect of plant density on plant growth in height.
        # Call add_plant_height to compute PlantHeight.
        if self.version < 0x500 or not state.date >= self.topping_date:
            # node numbers of top node.

            if len(state.vegetative_branches[0].fruiting_branches) >= 2:
                stage = state.vegetative_branches[0].fruiting_branches[1].nodes[0].stage
            else:
                stage = Stage.NotYetFormed
            state.plant_height += state.add_plant_height(
                denf2,
                state.day_inc,
                state.number_of_pre_fruiting_nodes,
                stage,
                state.age_of_pre_fruiting_nodes[state.number_of_pre_fruiting_nodes - 1],
                state.age_of_pre_fruiting_nodes[state.number_of_pre_fruiting_nodes - 2],
                state.agetop,
                state.water_stress_stem,
                state.carbon_stress,
                state.nitrogen_stress_vegetative,
                *self.cultivar_parameters[19:27]
            )
        # Call ActualRootGrowth() to compute actual root growth.
        state.compute_actual_root_growth(
            sumpdr,
            self.row_space,
            self.per_plant_area,
            3,
            self.emerge_date,
            self.plant_row_column,
        )

    @property
    def _column_width(self):
        return self.row_space / 20
