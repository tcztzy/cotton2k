cdef extern from "Simulation.h":
    struct Simulation:
        pass

cdef extern from "Cottonmodel.h":

    cdef cppclass C2KApp:
        C2KApp() except +
        void RunTheModel(const char *)
        void DailySimulation(Simulation &)
        void SimulateThisDay(Simulation &, const int &)

def run(str profile):
    app = new C2KApp()
    app.RunTheModel(profile.encode("utf-8"))
