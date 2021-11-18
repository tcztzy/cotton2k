cmake_minimum_required(VERSION 3.17)

project(cotton2k)
# Use CXX 17 has beta built-in function.
set(CMAKE_CXX_STANDARD 14)
file(GLOB SOURCE_FILES "*.cpp" "*.h")
add_library(cotton2k STATIC ${SOURCE_FILES})