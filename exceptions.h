#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H
#include <stdexcept>

class SimulationEnd
    : public std::runtime_error
{
public:
    SimulationEnd() noexcept
        : std::runtime_error("Simulation end.")
    {
    }

    SimulationEnd(const std::string &_Message) noexcept
        : std::runtime_error(_Message)
    {
    }
};
#endif //EXCEPTIONS_H