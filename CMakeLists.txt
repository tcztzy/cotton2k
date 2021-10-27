cmake_minimum_required(VERSION 3.17)

project(cotton2k)
# Use CXX 17 has beta built-in function.
set(CMAKE_CXX_STANDARD 14)
find_package(Boost 1.77.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
file(GLOB SOURCE_FILES "*.cpp" "*.h")
add_library(cotton2k STATIC ${SOURCE_FILES})