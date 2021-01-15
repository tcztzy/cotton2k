cmake_minimum_required(VERSION 3.17)
project(cotton2k)
# Use CXX 17 for filesystem and so on.
set(CMAKE_CXX_STANDARD 17)
# Use MFC in a Static Library
add_definitions(-D_AFXDLL)
set(CMAKE_MFC_FLAG 2)
file(GLOB SOURCE_FILES "*.cpp" "*.h" "Cottonmodel.rc")
set_source_files_properties(Cottonmodel.rc PROPERTIES LANGUAGE RC)
add_executable(cotton2k WIN32 ${SOURCE_FILES})
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
add_dependencies(cotton2k cotton2k.rs)
target_link_libraries(cotton2k
    debug "${CMAKE_SOURCE_DIR}/target/debug/cotton2k.lib"
    optimized "${CMAKE_SOURCE_DIR}/target/release/cotton2k.lib"
    ws2_32 userenv advapi32)
