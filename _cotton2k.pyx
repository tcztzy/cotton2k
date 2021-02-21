cdef extern from "Simulation.h":
    struct Simulation:
        pass

cdef extern from "Input.h":
    Simulation ReadInput(const char *)

cdef extern from "Output.h":
    void DataOutput(Simulation &)

cdef extern:
    void output_json(const Simulation &)

cdef extern from "Cottonmodel.h":

    cdef cppclass C2KApp:
        C2KApp() except +
        void RunTheModel(const char *)
        void DailySimulation(Simulation &)
        void SimulateThisDay(Simulation &, const int &)

def run(str profile):
    app = new C2KApp()
    sim = ReadInput(profile.encode("utf-8"))
    app.DailySimulation(sim);
    DataOutput(sim)
    output_json(sim)
