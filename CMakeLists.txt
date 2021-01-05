cmake_minimum_required(VERSION 3.17)
project(cotton2k)
# Use CXX 17 for filesystem and so on.
set(CMAKE_CXX_STANDARD 17)
# Use MFC in a Static Library
add_definitions(-D_AFXDLL)
set(CMAKE_MFC_FLAG 2)
set_source_files_properties(Cottonmodel.rc PROPERTIES LANGUAGE RC)
add_executable(cotton2k WIN32
        Cottonmodel.h Cottonmodel.cpp Cottonmodel.rc
        Dialogs.h Dialogs.cpp
        MainFrm.h MainFrm.cpp
        GeneralFunctions.h GeneralFunctions.cpp
        CottonPhenology.h CottonPhenology.cpp
        DailyClimate.h DailyClimate.cpp
        Input.h GettingInput_1.cpp GettingInput_2.cpp gettingInput_3.cpp
        Output.h WriteOutput_1.cpp WriteOutput_2.cpp WriteOutput_3.cpp
        PlantAdjustment.h PlantAdjustment.cpp
        PlantGrowth.h PlantGrowth_1.cpp PlantGrowth_2.cpp PlantGrowth_3.cpp
        PlantNitrogen.h PlantNitrogen.cpp
        SoilNitrogen.h SoilNitrogen.cpp
        SoilProcedures.h SoilProcedures_1.cpp SoilProcedures_2.cpp SoilProcedures_3.cpp
        SoilTemperature.h SoilTemperature_1.cpp SoilTemperature_2.cpp SoilTemperature_3.cpp
        global.h global.cpp InitializeGlobal.cpp
        LeafAbscission.h LeafAbscission.cpp
        FruitAbscission.h FruitAbscission.cpp
        RootGrowth.h RootGrowth_1.cpp RootGrowth_2.cpp)
# Enable ExternalProject CMake module
include(ExternalProject)
ExternalProject_Add(
        cotton2k.rs
        DOWNLOAD_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND cargo build COMMAND cargo build --release
        BINARY_DIR "${CMAKE_SOURCE_DIR}"
        INSTALL_COMMAND ""
        LOG_BUILD ON
)
add_dependencies(cotton2k cotton2k.rs)
target_link_libraries(cotton2k
    debug "${CMAKE_SOURCE_DIR}/target/debug/cotton2k.lib"
    optimized "${CMAKE_SOURCE_DIR}/target/release/cotton2k.lib"
    ws2_32 userenv advapi32)
