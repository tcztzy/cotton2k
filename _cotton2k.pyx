from functools import cache
from datetime import date, timedelta

cdef extern from "State.h":
    ctypedef struct State:
        double plant_height
        double plant_weight
        double lint_yield

cdef extern from "Simulation.h":
    ctypedef struct Simulation:
        int year
        unsigned int day_start
        unsigned int day_finish
        unsigned int day_emerge
        State *states

cdef extern from "Input.h":
    Simulation ReadInput(const char *)

cdef extern from "Output.h":
    void DataOutput(Simulation &)

cdef extern from "Cottonmodel.h":

    cdef cppclass C2KApp:
        C2KApp() except +
        void RunTheModel(const char *)
        void DailySimulation(Simulation &)
        void SimulateThisDay(Simulation &, const int &)

cdef class _Simulation:
    cdef Simulation _sim

    @property
    def start_date(self):
        return date(self._sim.year, 1, 1) + timedelta(days=self._sim.day_start - 1)

    @property
    def end_date(self):
        return date(self._sim.year, 1, 1) + timedelta(days=self._sim.day_finish - 1)

    @property
    def emerge_date(self):
        return date(self._sim.year, 1, 1) + timedelta(days=self._sim.day_emerge - 1)

    @property
    def states(self):
        return (self._sim.states[i] for i in range(self._sim.day_finish - self._sim.day_start + 1))

def run(str profile):
    app = new C2KApp()
    sim = ReadInput(profile.encode("utf-8"))
    app.DailySimulation(sim);
    DataOutput(sim)
    py_sim = _Simulation()
    py_sim._sim = sim
    return py_sim
