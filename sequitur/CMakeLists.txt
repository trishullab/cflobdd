project(sequitur)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wfatal-errors -std=c++11 -flto")

set(PROJECT_SOURCES
    ${CMAKE_SOURCE_DIR}/main.cpp
    ${CMAKE_SOURCE_DIR}/sequitur.hpp

    ${CMAKE_SOURCE_DIR}/sequitur/baselist.hpp
    ${CMAKE_SOURCE_DIR}/sequitur/baselist.tpp
    ${CMAKE_SOURCE_DIR}/sequitur/hashing.hpp
    ${CMAKE_SOURCE_DIR}/sequitur/id.hpp
    ${CMAKE_SOURCE_DIR}/sequitur/objectpool.hpp
    ${CMAKE_SOURCE_DIR}/sequitur/symbols.cpp
    ${CMAKE_SOURCE_DIR}/sequitur/symbols.hpp
    ${CMAKE_SOURCE_DIR}/sequitur/symbolwrapper.cpp
    ${CMAKE_SOURCE_DIR}/sequitur/symbolwrapper.hpp
)

include_directories(
    ${CMAKE_SOURCE_DIR} sequitur
)

add_executable(${PROJECT_NAME} ${PROJECT_SOURCES})
