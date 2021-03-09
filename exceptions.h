#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H
#include <stdexcept>
#include <filesystem>
namespace fs = std::filesystem;

class Cotton2KException
    : public std::runtime_error
{
public:
    Cotton2KException(const std::string &_Message)
        : std::runtime_error(_Message.c_str())
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
    : public std::runtime_error
{
public:
    SimulationEnd() noexcept
        : std::runtime_error("Simulation end.")
    {
    }

    SimulationEnd(const std::string &_Message) noexcept
        : std::runtime_error(_Message.c_str())
    {
    }
};
#endif //EXCEPTIONS_H