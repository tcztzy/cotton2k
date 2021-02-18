cmake_minimum_required(VERSION 3.17)
project(cotton2k VERSION 5.0.0 LANGUAGES CXX C)

# Use CXX 17 for filesystem and so on.
set(CMAKE_CXX_STANDARD 17)
# Use MFC in a Static Library
add_definitions(-D_AFXDLL)
set(CMAKE_MFC_FLAG 2)

# Define the SOURCE_FILES
file(GLOB SOURCE_FILES "*.cpp" "*.h" "cotton2k.pyx")

# Prefer virtualenv
set(Python_FIND_VIRTUALENV ONLY)
# Sci-kit Build
find_package(PythonExtensions REQUIRED)
find_package(Cython REQUIRED)

set(LIBRARY _${PROJECT_NAME})

add_cython_target(${LIBRARY} CXX)
add_library(${LIBRARY} MODULE ${${LIBRARY}} ${SOURCE_FILES})
python_extension_module(${LIBRARY})

target_compile_options(${LIBRARY} PRIVATE
        "$<$<CXX_COMPILER_ID:Clang>:-Weverything;-Wfatal-errors>"
        "$<$<CXX_COMPILER_ID:GNU>:-Wall;-Wextra;-Wfatal-errors>"
        "$<$<CXX_COMPILER_ID:MSVC>:/W3>")

# Enable ExternalProject CMake module
include(ExternalProject)
ExternalProject_Add(cotton2k.rs
        DOWNLOAD_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND cargo build COMMAND cargo build --release
        BUILD_ALWAYS 1
        BINARY_DIR "${CMAKE_SOURCE_DIR}"
        INSTALL_COMMAND ""
        LOG_BUILD ON
)
add_dependencies(${LIBRARY} cotton2k.rs)
target_link_libraries(${LIBRARY}
        debug "${CMAKE_SOURCE_DIR}/target/debug/cotton2k.lib"
        optimized "${CMAKE_SOURCE_DIR}/target/release/cotton2k.lib"
        ws2_32 userenv advapi32)
target_include_directories(${LIBRARY} PRIVATE
        ${CMAKE_SOURCE_DIR})
install(TARGETS ${LIBRARY} LIBRARY DESTINATION cotton2k)
