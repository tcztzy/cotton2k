# distutils: language=c++
# cython: language_level=3
import datetime

import numpy as np

from _cotton2k.utils import date2doy, doy2date

cdef class StateBase:
    rlat1 = np.zeros(40, dtype=np.float64)  # lateral root length (cm) to the left of the tap root
    rlat2 = np.zeros(40, dtype=np.float64)   # lateral root length (cm) to the right of the tap root
    actual_root_growth = np.zeros((40, 20), dtype=np.float64)
    root_potential_growth = np.zeros((40, 20), dtype=np.float64)  # potential root growth in a soil cell (g per day).
    root_age = np.zeros((40, 20), dtype=np.float64)

    def keys(self):
        return [
            "lint_yield",
            "ginning_percent",
            "leaf_area_index",
            "light_interception",
            "plant_height",
            *(org + "_weight" for org in ("plant", "stem", "leaf", "root", "petiole", "square", "green_bolls", "open_bolls")),
            "evapotranspiration",
            "actual_transpiration",
            "actual_soil_evaporation",
            *("number_of_" + org for org in ("squares", "green_bolls", "open_bolls")),
            "vegetative_branches",
        ]

    @property
    def daynum(self):
        return self._[0].daynum

    @daynum.setter
    def daynum(self, value):
        self._[0].daynum = value

    @property
    def date(self):
        return doy2date(self.year, self._[0].daynum)

    @date.setter
    def date(self, value):
        if not isinstance(value, datetime.date):
            raise TypeError
        self._[0].daynum = date2doy(value)

    @property
    def day_inc(self):
        return self._[0].day_inc

    @day_inc.setter
    def day_inc(self, value):
        self._[0].day_inc = value

    @property
    def kday(self):
        return self._[0].kday

    @kday.setter
    def kday(self, value):
        self._[0].kday = value

    @property
    def average_temperature(self):
        return self._[0].average_temperature

    @average_temperature.setter
    def average_temperature(self, value):
        self._[0].average_temperature = value

    @property
    def lint_yield(self):
        return self._[0].lint_yield

    @lint_yield.setter
    def lint_yield(self, value):
        self._[0].lint_yield = value

    @property
    def plant_height(self):
        return self._[0].plant_height

    @plant_height.setter
    def plant_height(self, value):
        self._[0].plant_height = value

    @property
    def plant_weight(self):
        return self._[0].plant_weight

    @plant_weight.setter
    def plant_weight(self, value):
        self._[0].plant_weight = value

    @property
    def stem_weight(self):
        return self._[0].stem_weight

    @stem_weight.setter
    def stem_weight(self, value):
        self._[0].stem_weight = value

    @property
    def petiole_weight(self):
        return self._[0].petiole_weight

    @petiole_weight.setter
    def petiole_weight(self, value):
        self._[0].petiole_weight = value

    @property
    def square_weight(self):
        return self._[0].square_weight

    @square_weight.setter
    def square_weight(self, value):
        self._[0].square_weight = value

    @property
    def green_bolls_weight(self):
        return self._[0].green_bolls_weight

    @green_bolls_weight.setter
    def green_bolls_weight(self, value):
        self._[0].green_bolls_weight = value

    @property
    def green_bolls_burr_weight(self):
        return self._[0].green_bolls_burr_weight

    @green_bolls_burr_weight.setter
    def green_bolls_burr_weight(self, value):
        self._[0].green_bolls_burr_weight = value

    @property
    def open_bolls_weight(self):
        return self._[0].open_bolls_weight

    @open_bolls_weight.setter
    def open_bolls_weight(self, value):
        self._[0].open_bolls_weight = value

    @property
    def open_bolls_burr_weight(self):
        return self._[0].open_bolls_burr_weight

    @open_bolls_burr_weight.setter
    def open_bolls_burr_weight(self, value):
        self._[0].open_bolls_burr_weight = value

    @property
    def root_weight(self):
        return self._[0].root_weight

    @root_weight.setter
    def root_weight(self, value):
        self._[0].root_weight = value

    @property
    def bloom_weight_loss(self):
        return self._[0].bloom_weight_loss

    @bloom_weight_loss.setter
    def bloom_weight_loss(self, value):
        self._[0].bloom_weight_loss = value

    @property
    def abscised_fruit_sites(self):
        return self._[0].abscised_fruit_sites

    @abscised_fruit_sites.setter
    def abscised_fruit_sites(self, value):
        self._[0].abscised_fruit_sites = value

    @property
    def abscised_leaf_weight(self):
        return self._[0].abscised_leaf_weight

    @abscised_leaf_weight.setter
    def abscised_leaf_weight(self, value):
        self._[0].abscised_leaf_weight = value

    @property
    def cumulative_nitrogen_loss(self):
        return self._[0].cumulative_nitrogen_loss

    @cumulative_nitrogen_loss.setter
    def cumulative_nitrogen_loss(self, value):
        self._[0].cumulative_nitrogen_loss = value

    @property
    def cumulative_transpiration(self):
        return self._[0].cumulative_transpiration

    @cumulative_transpiration.setter
    def cumulative_transpiration(self, value):
        self._[0].cumulative_transpiration = value

    @property
    def cumulative_evaporation(self):
        return self._[0].cumulative_evaporation

    @cumulative_evaporation.setter
    def cumulative_evaporation(self, value):
        self._[0].cumulative_evaporation = value

    @property
    def applied_water(self):
        return self._[0].applied_water

    @applied_water.setter
    def applied_water(self, value):
        self._[0].applied_water = value

    @property
    def water_stress(self):
        return self._[0].water_stress

    @water_stress.setter
    def water_stress(self, value):
        self._[0].water_stress = value

    @property
    def water_stress_stem(self):
        return self._[0].water_stress_stem

    @water_stress_stem.setter
    def water_stress_stem(self, value):
        self._[0].water_stress_stem = value

    @property
    def carbon_stress(self):
        return self._[0].carbon_stress

    @carbon_stress.setter
    def carbon_stress(self, value):
        self._[0].carbon_stress = value

    @property
    def extra_carbon(self):
        return self._[0].extra_carbon

    @extra_carbon.setter
    def extra_carbon(self, value):
        self._[0].extra_carbon = value

    @property
    def light_interception(self):
        return self._[0].light_interception

    @light_interception.setter
    def light_interception(self, value):
        self._[0].light_interception = value

    @property
    def leaf_area_index(self):
        return self._[0].leaf_area_index

    @leaf_area_index.setter
    def leaf_area_index(self, value):
        self._[0].leaf_area_index = value

    @property
    def leaf_area(self):
        return self._[0].leaf_area

    @leaf_area.setter
    def leaf_area(self, value):
        self._[0].leaf_area = value

    @property
    def leaf_weight(self):
        return self._[0].leaf_weight

    @leaf_weight.setter
    def leaf_weight(self, value):
        self._[0].leaf_weight = value

    @property
    def leaf_weight_area_ratio(self):
        return self._[0].leaf_weight_area_ratio

    @leaf_weight_area_ratio.setter
    def leaf_weight_area_ratio(self, value):
        self._[0].leaf_weight_area_ratio = value

    @property
    def leaf_nitrogen(self):
        return self._[0].leaf_nitrogen

    @leaf_nitrogen.setter
    def leaf_nitrogen(self, value):
        self._[0].leaf_nitrogen = value

    @property
    def leaf_nitrogen_concentration(self):
        return self._[0].leaf_nitrogen_concentration

    @leaf_nitrogen_concentration.setter
    def leaf_nitrogen_concentration(self, value):
        self._[0].leaf_nitrogen_concentration = value

    @property
    def number_of_vegetative_branches(self):
        return self._[0].number_of_vegetative_branches

    @number_of_vegetative_branches.setter
    def number_of_vegetative_branches(self, value):
        self._[0].number_of_vegetative_branches = value

    @property
    def number_of_squares(self):
        return self._[0].number_of_squares

    @number_of_squares.setter
    def number_of_squares(self, value):
        self._[0].number_of_squares = value

    @property
    def number_of_green_bolls(self):
        return self._[0].number_of_green_bolls

    @number_of_green_bolls.setter
    def number_of_green_bolls(self, value):
        self._[0].number_of_green_bolls = value

    @property
    def number_of_open_bolls(self):
        return self._[0].number_of_open_bolls

    @number_of_open_bolls.setter
    def number_of_open_bolls(self, value):
        self._[0].number_of_open_bolls = value

    @property
    def nitrogen_stress(self):
        return self._[0].nitrogen_stress

    @nitrogen_stress.setter
    def nitrogen_stress(self, value):
        self._[0].nitrogen_stress = value

    @property
    def nitrogen_stress_vegetative(self):
        return self._[0].nitrogen_stress_vegetative

    @nitrogen_stress_vegetative.setter
    def nitrogen_stress_vegetative(self, value):
        self._[0].nitrogen_stress_vegetative = value

    @property
    def nitrogen_stress_fruiting(self):
        return self._[0].nitrogen_stress_fruiting

    @nitrogen_stress_fruiting.setter
    def nitrogen_stress_fruiting(self, value):
        self._[0].nitrogen_stress_fruiting = value

    @property
    def nitrogen_stress_root(self):
        return self._[0].nitrogen_stress_root

    @nitrogen_stress_root.setter
    def nitrogen_stress_root(self, value):
        self._[0].nitrogen_stress_root = value

    @property
    def total_required_nitrogen(self):
        return self._[0].total_required_nitrogen

    @total_required_nitrogen.setter
    def total_required_nitrogen(self, value):
        self._[0].total_required_nitrogen = value

    @property
    def petiole_nitrogen_concentration(self):
        return self._[0].petiole_nitrogen_concentration

    @petiole_nitrogen_concentration.setter
    def petiole_nitrogen_concentration(self, value):
        self._[0].petiole_nitrogen_concentration = value

    @property
    def seed_nitrogen(self):
        return self._[0].seed_nitrogen

    @seed_nitrogen.setter
    def seed_nitrogen(self, value):
        self._[0].seed_nitrogen = value

    @property
    def seed_nitrogen_concentration(self):
        return self._[0].seed_nitrogen_concentration

    @seed_nitrogen_concentration.setter
    def seed_nitrogen_concentration(self, value):
        self._[0].seed_nitrogen_concentration = value

    @property
    def burr_nitrogen(self):
        return self._[0].burr_nitrogen

    @burr_nitrogen.setter
    def burr_nitrogen(self, value):
        self._[0].burr_nitrogen = value

    @property
    def root_nitrogen_concentration(self):
        return self._[0].root_nitrogen_concentration

    @root_nitrogen_concentration.setter
    def root_nitrogen_concentration(self, value):
        self._[0].root_nitrogen_concentration = value

    @property
    def square_nitrogen(self):
        return self._[0].square_nitrogen

    @square_nitrogen.setter
    def square_nitrogen(self, value):
        self._[0].square_nitrogen = value

    @property
    def square_nitrogen_concentration(self):
        return self._[0].square_nitrogen_concentration

    @square_nitrogen_concentration.setter
    def square_nitrogen_concentration(self, value):
        self._[0].square_nitrogen_concentration = value

    @property
    def stem_nitrogen(self):
        return self._[0].stem_nitrogen

    @stem_nitrogen.setter
    def stem_nitrogen(self, value):
        self._[0].stem_nitrogen = value

    @property
    def fruit_growth_ratio(self):
        return self._[0].fruit_growth_ratio

    @fruit_growth_ratio.setter
    def fruit_growth_ratio(self, value):
        self._[0].fruit_growth_ratio = value

    @property
    def ginning_percent(self):
        return self._[0].ginning_percent

    @ginning_percent.setter
    def ginning_percent(self, value):
        self._[0].ginning_percent = value

    @property
    def number_of_pre_fruiting_nodes(self):
        return self._[0].number_of_pre_fruiting_nodes

    @number_of_pre_fruiting_nodes.setter
    def number_of_pre_fruiting_nodes(self, value):
        self._[0].number_of_pre_fruiting_nodes = value

    @property
    def age_of_pre_fruiting_nodes(self):
        return self._[0].age_of_pre_fruiting_nodes

    @property
    def leaf_area_pre_fruiting(self):
        return self._[0].leaf_area_pre_fruiting

    @property
    def leaf_weight_pre_fruiting(self):
        return self._[0].leaf_weight_pre_fruiting

    @property
    def evapotranspiration(self):
        return self._[0].evapotranspiration

    @evapotranspiration.setter
    def evapotranspiration(self, value):
        self._[0].evapotranspiration = value

    @property
    def actual_transpiration(self):
        return self._[0].actual_transpiration

    @actual_transpiration.setter
    def actual_transpiration(self, value):
        self._[0].actual_transpiration = value

    @property
    def actual_soil_evaporation(self):
        return self._[0].actual_soil_evaporation

    @actual_soil_evaporation.setter
    def actual_soil_evaporation(self, value):
        self._[0].actual_soil_evaporation = value

    @property
    def pavail(self):
        return self._[0].pavail

    @pavail.setter
    def pavail(self, value):
        self._[0].pavail = value

    @property
    def deep_soil_temperature(self):
        return self._[0].deep_soil_temperature

    @deep_soil_temperature.setter
    def deep_soil_temperature(self, value):
        self._[0].deep_soil_temperature = value

    def __getitem__(self, item):
        return getattr(self, item)
