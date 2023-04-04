cmake_minimum_required(VERSION 3.16)
project(ska_cpp)
set(CMAKE_CXX_STANDARD 17)

# Variable definitions
set(TARGET_NAME ska_cpp)
add_compile_definitions(PYTHON_EXT)

set(CMAKE_OSX_ARCHITECTURES "x86_64")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=apple-m1")

# gcc: Add openmp
#IF(CMAKE_COMPILER_IS_GNUCC OR "$ENV{SKETCHLIB_INSTALL}" STREQUAL "conda")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
#ENDIF()
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g")
# Add -O0 to remove optimizations when using debug
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

if(UNIX AND NOT APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS")
    set(CMAKE_LD_FLAGS "${CMAKE_LDFLAGS} -Wl,--as-needed")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math -funroll-loops -m64")

# Set paths for non standard lib/ and include/ locations
if(DEFINED ENV{CONDA_PREFIX})
    include_directories($ENV{CONDA_PREFIX}/include)
    link_directories($ENV{CONDA_PREFIX}/lib)
    link_directories($ENV{CONDA_PREFIX}/lib/intel64)
endif()

# Add libraries

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/cereal/include)


find_package(ZLIB)
execute_process(COMMAND pybind11-config --cmakedir OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE pybind11_DIR)
find_package(pybind11 2.6 CONFIG REQUIRED)

# Define python library target
add_library("${TARGET_NAME}" MODULE src/old_code.cpp)

# Compile CPU library
target_sources("${TARGET_NAME}" PRIVATE
    src/ska_fasta.cpp
    src/ska_align.cpp
    src/ska_bindings.cpp)

set_target_properties("${TARGET_NAME}" PROPERTIES
    CXX_VISIBILITY_PRESET "hidden"
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
)
if(UNIX AND (NOT APPLE OR NOT CMAKE_COMPILER_IS_GNUCC))
    set_target_properties("${TARGET_NAME}" PROPERTIES
                          INTERPROCEDURAL_OPTIMIZATION ON)
endif()

# Link libraries
target_link_libraries("${TARGET_NAME}" PRIVATE pybind11::module m dl)
if(DEFINED ENV{CONDA_PREFIX} AND (NOT APPLE OR CMAKE_COMPILER_IS_GNUCC OR ENV{SKETCHLIB_INSTALL} EQUAL "conda"))
    target_link_libraries("${TARGET_NAME}" PRIVATE gomp z)
else()
    target_link_libraries("${TARGET_NAME}" PRIVATE ZLIB::ZLIB)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries("${TARGET_NAME}" PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()
