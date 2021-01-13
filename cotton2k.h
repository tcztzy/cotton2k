#ifndef COTTON2K_H
#define COTTON2K_H
#include <filesystem>
namespace fs = std::filesystem;
typedef struct ClimateStruct {
    int nDay;
    double Rad, Tmax, Tmin, Rain, Wind, Tdew;
} ClimateStruct;

class Cotton2KException
        : public std::exception {
public:
    Cotton2KException(const std::string &_Message)
            : exception(_Message.c_str(), 1) {
    }
};

class FileNotExists
        : public Cotton2KException {
public:

    FileNotExists(const fs::path filePath) noexcept
            : Cotton2KException("The file  " + filePath.string() + "  does not exist!") {
    }
};

class FileNotOpened
        : public Cotton2KException {
public:

    FileNotOpened(const fs::path filePath) noexcept
            : Cotton2KException("Can't open " + filePath.string() + ".") {
    }
};

class SimulationEnd
        : public std::exception {
public:

    SimulationEnd() noexcept
            : std::exception("Simulation end.", 1) {
    }

    SimulationEnd(const std::string &_Message) noexcept
            : std::exception(_Message.c_str(), 1) {
    }
};
#endif//COTTON2K_H