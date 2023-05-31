cmake_minimum_required(VERSION 3.14)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Build type")

set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release)

# Organisation of Visual Studio solution
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY AUTOGEN_TARGETS_FOLDER AutoMoc)
set_property(GLOBAL PROPERTY AUTOGEN_SOURCE_GROUP AutoMoc)
set(CMAKE_UtilityFolderName "CMakeUtilityTargets")
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER ${CMAKE_UtilityFolderName})

set(IPP_ROOT CACHE PATH "IPP library path") 

if ( NOT IPP_ROOT )
message(FATAL_ERROR "You need to provide the correct path to all specified libraries and sources."
        "It is assumed that libraries are"
        "compartable with CMake, and `find_package(<Path>)` can correctly process configuration."     
        "`<Lib>Config.cmake` must be under that path."
        "Sources must have CMakeLists.txt"
)
endif()

project(KeysightTechnologiesControl C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(QT_ROOT "Path to QT_ROOT" CACHE PATH "Path to Qt installed location.")

find_path(QT_CMAKE_PATH Qt5Config.cmake ${QT_ROOT}/lib/cmake/Qt5)
mark_as_advanced(QT_CMAKE_PATH)

if(NOT QT_CMAKE_PATH)
    message(FATAL_ERROR 
        "QT installation not found. Please set QT_ROOT to point to the QT installation.
        e.g C:/Qt/Qt5.9.1/5.9.1/msvc2015_64"
    )
endif()

# Instruct CMake to run moc automatically when needed
set(CMAKE_AUTOMOC ON)
# Create code from a list of Qt designer ui files
set(CMAKE_AUTOUIC ON)
# Create code to handle resources
set(CMAKE_AUTORCC ON)

list(APPEND CMAKE_PREFIX_PATH ${QT_ROOT})

find_package(Qt5 
    COMPONENTS 
        Widgets 
        DataVisualization
        Core
        LinguistTools
        Gui
        Network
        PrintSupport
        Concurrent
        Svg
    REQUIRED
)

add_executable(KeysightTechnologiesControl src/main.cpp)

if(MSVC)
    target_compile_options(KeysightTechnologiesControl PRIVATE "/MP")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(IPP REQUIRED)

target_compile_definitions(KeysightTechnologiesControl PRIVATE 
    _USE_MATH_DEFINES
    _UNICODE
    _NOMINMAX
    NO_OPTIX
)

target_sources(KeysightTechnologiesControl PRIVATE
    src/main.cpp
    src/sentinel1_packet_decode.h
    src/sentinel1_packet_decode.cpp
    src/sentinel1.h
    src/sentinel1.cpp
    src/RLI.cpp
    src/RLI.h
    src/RLI.ui
    src/zoom.cpp
    src/zoom.h
    src/radarsat1.h
    src/radarsat1.cpp
    src/qcustomplot.cpp
    src/qcustomplot.h
    src/volume_dir_file.h
    src/trl_file.h
    src/dat_file.h
    src/null_vol_dir_file.h
    src/ldr_file.h
)

target_include_directories(KeysightTechnologiesControl PRIVATE  
    ${IPP_INCLUDE_DIR}
)

target_link_libraries(KeysightTechnologiesControl PRIVATE
    ipp::core
    ipp::i
    ipp::s
    ipp::vm
)

add_custom_command(
    TARGET KeysightTechnologiesControl POST_BUILD  
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${IPP_DLLS} $<TARGET_FILE_DIR:KeysightTechnologiesControl>
) 

target_link_libraries(KeysightTechnologiesControl PRIVATE
    Qt5::Concurrent
    Qt5::Core
    Qt5::DataVisualization
    Qt5::Gui
    Qt5::Network
    Qt5::PrintSupport
    Qt5::Widgets
    Qt5::Svg
    #Boost::log
    #Boost::system
    #Boost::log_setup
)

target_include_directories(KeysightTechnologiesControl PRIVATE "Eigen")

add_custom_command(
    TARGET KeysightTechnologiesControl POST_BUILD  
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_FILE:Qt5::Core>
        $<TARGET_FILE:Qt5::DataVisualization>
        $<TARGET_FILE:Qt5::Gui>
        $<TARGET_FILE:Qt5::Network>
        $<TARGET_FILE:Qt5::PrintSupport>
        $<TARGET_FILE:Qt5::Widgets>
        $<TARGET_FILE:Qt5::Svg>
        $<TARGET_FILE_DIR:KeysightTechnologiesControl>
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${QT_ROOT}/plugins/platforms
        $<TARGET_FILE_DIR:KeysightTechnologiesControl>/platforms
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${IPP_DLLS} $<TARGET_FILE_DIR:KeysightTechnologiesControl>
) 