cdef extern from "Irrigation.h":
    ctypedef struct Irrigation:
        int day
        int method
        int LocationColumnDrip
        int LocationLayerDrip
        double amount
