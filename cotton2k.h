#ifndef COTTON2K_H
#define COTTON2K_H
#include <filesystem>
namespace fs = std::filesystem;
typedef struct ClimateStruct
{
    int nDay;
    double Rad, Tmax, Tmin, Rain, Wind, Tdew;
} ClimateStruct;

typedef struct State
{
    uint32_t day_number = 0; // Days from the start of the first year of simulation (day of year = DOY)
    char date[12];           // date string formatted as "dd-MMM-yyyy", for example 25-JUN-2003
    double plant_height = 4.0;     // plant height, cm.
} State;

typedef struct Simulation
{
    const char *profile_name;               // name of input file with profile data (without the extension ".PRO")
    size_t profile_name_length;             //
    int day_emerge;                         // Date of emergence (DOY).
    int day_start;                          // Date (DOY) to start simulation.
    int day_finish;                         // Date (DOY) to finish simulation.
    int day_plant;                          // Date (DOY) of planting.
    int day_start_soil_maps;                // Date (DOY) to start soil slab maps output.
    int day_stop_soil_maps;                 // Date (DOY) to stop soil slab maps output.
    int day_start_co2;                      // First date (DOY) with CO2 enrichment.
    int day_end_co2;                        // Last date (DOY) with CO2 enrichment.
    double co2_enrichment_factor;           // factor describing effect of CO2 enrichment.
    int day_start_mulch;                    // Date (DOY) for beginning of mulch.
    int day_end_mulch;                      // date (DOY) for ending of mulch.
    int mulch_indicator;                    // indicating if and where a soil mulch exists, the value are:
                                            // 0 = no mulch;
                                            // 1 = plastic layer on all soil surface;
                                            // 2 = plastic layer on all soil surface except one column at each side of the plant row;
                                            // 3 = plastic layer on all soil surface except two columns at each side of the plant row.
    double mulch_transmissivity_short_wave; // transmissivity of soil mulch to short wave radiation
    double mulch_transmissivity_long_wave;  // transmissivity of soil mulch to long wave radiation.
    double latitude;                        // degree
    double longitude;                       // degree
    int num_curve;                          // number of input soil-moisture curves in the impedance table.
    int first_bloom = 0;                    // Date (DOY) of first bloom.
    int first_square = 0;                   // Date of first square (DOY), if no squares have been formed, FirstSquare = 0.
    double root_weight[40][20][3];          // weight of dry matter of root tissue in a soil cell for an age group, in g per cell.
    double root_age[40][20];                // the time (in days) from the first appearance of roots in a soil cell.
    ClimateStruct climate[400];             // structure containing the following daily weather data:
                                            // int nDay =    day of year.
                                            // double Rad =  daily global radiation, in langleys.
                                            // double Tmax = maximum daily temperature, C.
                                            // double Tmin = minimum daily temperature, C.
                                            // double Tdew = dew point temperature, C.
                                            // double Rain = daily rainfall, mm.
                                            // double Wind = daily wind run, km.
    State *states;
    size_t number_of_states;
} Simulation;

class Cotton2KException
    : public std::exception
{
public:
    Cotton2KException(const std::string &_Message)
        : exception(_Message.c_str(), 1)
    {
    }
};

class FileNotExists
    : public Cotton2KException
{
public:
    FileNotExists(const fs::path filePath) noexcept
        : Cotton2KException("The file  " + filePath.string() + "  does not exist!")
    {
    }
};

class FileNotOpened
    : public Cotton2KException
{
public:
    FileNotOpened(const fs::path filePath) noexcept
        : Cotton2KException("Can't open " + filePath.string() + ".")
    {
    }
};

class SimulationEnd
    : public std::exception
{
public:
    SimulationEnd() noexcept
        : std::exception("Simulation end.", 1)
    {
    }

    SimulationEnd(const std::string &_Message) noexcept
        : std::exception(_Message.c_str(), 1)
    {
    }
};
#endif //COTTON2K_H