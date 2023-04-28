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

add_executable(KeysightTechnologiesControl src/main.cpp)
 
if(MSVC)
    target_compile_options(KeysightTechnologiesControl PRIVATE "/MP")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(IPP REQUIRED)

add_subdirectory(boost)

target_compile_definitions(KeysightTechnologiesControl PRIVATE 
    _USE_MATH_DEFINES
    _UNICODE
    _NOMINMAX
    NO_OPTIX
)

target_sources(KeysightTechnologiesControl PRIVATE
    src/main.cpp
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
    Boost::log
    Boost::system
    Boost::log_setup
)

include(FetchContent)
    
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
    
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
    
target_link_libraries(KeysightTechnologiesControl PRIVATE GTest::gtest_main)
    
include(GoogleTest)

