cmake_minimum_required(VERSION 3.12)
project(cflobdd)

# Set the C++ standard to 17, as specified in the compile command
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add compiler flags: -g for debugging and -w to suppress warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -w")

# Define the main source directory
set(CFLOBDD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CFLOBDD)

# Add include directories.
# The user can specify the Boost path via CMAKE_PREFIX_PATH or by setting BOOST_ROOT.
find_package(Boost REQUIRED CONFIG)
include_directories(
    ${CFLOBDD_DIR}
    ${Boost_INCLUDE_DIRS}
    ${CFLOBDD_DIR}/Solver/uwr/bit_vector
    ${CFLOBDD_DIR}/Solver/uwr/assert
    ${CFLOBDD_DIR}/Solver/uwr/matrix
    ${CFLOBDD_DIR}/Solver/uwr/parsing
)

# Gather the source files.
file(GLOB SOURCES
    "${CFLOBDD_DIR}/*.cpp"
    "${CFLOBDD_DIR}/Solver/uwr/bit_vector/*.cpp"
    "${CFLOBDD_DIR}/Solver/uwr/parsing/*.cpp"
    "${CFLOBDD_DIR}/cflobdd_c.cpp"
)
list(REMOVE_ITEM SOURCES "${CFLOBDD_DIR}/main.cpp")

# Create the shared library.
add_library(cflobdd SHARED ${SOURCES})
add_executable(cflobdd_test ${SOURCES} ${CFLOBDD_DIR}/main.cpp)

# Link the math library.
target_link_libraries(cflobdd m)
target_compile_definitions(cflobdd PUBLIC PATH_COUNTING_ENABLED)

target_link_libraries(cflobdd_test m)
target_compile_definitions(cflobdd_test PUBLIC PATH_COUNTING_ENABLED)
