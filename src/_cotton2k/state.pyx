# distutils: language=c++
# cython: language_level=3
cdef class FruitingBranch:
    cdef cFruitingBranch _branch
    __slots__ = ("delay_for_new_node", "main_stem_leaf", "nodes")

    def __init__(self, _branch):
        self._branch = _branch

    @property
    def delay_for_new_node(self):
        return self._branch.delay_for_new_node

    @property
    def main_stem_leaf(self):
        return self._branch.main_stem_leaf

    @property
    def nodes(self):
        return [self._branch.nodes[i] for i in range(self._branch.number_of_fruiting_nodes)]

    def __iter__(self):
        for attr in self.__slots__:
            yield attr, getattr(self, attr)


cdef class VegetativeBranch:
    cdef cVegetativeBranch _branch

    def __init__(self, _branch):
        self._branch = _branch

    @property
    def fruiting_branches(self):
        return [FruitingBranch(self._branch.fruiting_branches[i]) for i in range(self._branch.number_of_fruiting_branches)]

    def __iter__(self):
        yield "fruiting_branches", self.fruiting_branches


cdef class State:
    cdef cState _state
    __slots__ = (
        "plant_height",
        "plant_weight",
        "lint_yield",
        "number_of_squares",
        "number_of_green_bolls",
        "number_of_open_bolls",
        "leaf_area_index",
        "ginning_percent",
        "vegetative_branches",
        "hours",
        "soil",
    )

    def __init__(self, _state):
        self._state = _state

    @property
    def daynum(self):
        return self._state.daynum

    @daynum.setter
    def daynum(self, value):
        self._state.daynum = value

    @property
    def plant_height(self):
        return self._state.plant_height

    @plant_height.setter
    def plant_height(self, value):
        self._state.plant_height = value

    @property
    def plant_weight(self):
        return self._state.plant_weight

    @plant_weight.setter
    def plant_weight(self, value):
        self._state.plant_weight = value

    @property
    def lint_yield(self):
        return self._state.lint_yield

    @lint_yield.setter
    def lint_yield(self, value):
        self._state.lint_yield = value

    @property
    def ginning_percent(self):
        return self._state.ginning_percent

    @ginning_percent.setter
    def ginning_percent(self, value):
        self._state.ginning_percent = value

    @property
    def number_of_squares(self):
        return self._state.number_of_squares

    @number_of_squares.setter
    def number_of_squares(self, value):
        self._state.number_of_squares = value

    @property
    def number_of_green_bolls(self):
        return self._state.number_of_green_bolls

    @number_of_green_bolls.setter
    def number_of_green_bolls(self, value):
        self._state.number_of_green_bolls = value

    @property
    def number_of_open_bolls(self):
        return self._state.number_of_open_bolls

    @number_of_open_bolls.setter
    def number_of_open_bolls(self, value):
        self._state.number_of_open_bolls = value

    @property
    def leaf_area_index(self):
        return self._state.leaf_area_index

    @leaf_area_index.setter
    def leaf_area_index(self, value):
        self._state.leaf_area_index = value

    @property
    def vegetative_branches(self):
        return [VegetativeBranch(self._state.vegetative_branches[i]) for i in range(self._state.number_of_vegetative_branches)]

    @property
    def hours(self):
        return self._state.hours

    @property
    def soil(self):
        return self._state.soil

    def __iter__(self):
        for attr in self.__slots__:
            value = getattr(self, attr)
            if value == 0 and attr.startswith("number_of_"):
                continue
            yield attr, value
