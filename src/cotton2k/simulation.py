# pylint: disable=no-name-in-module
from _cotton2k.simulation import Simulation as CySimulation


class Simulation(CySimulation):
    def read_input(self, *args, **kwargs):  # pylint: disable=unused-argument
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
            self._simulate_this_day(i)
            self._copy_state(i)
        self._simulate_this_day(days)
