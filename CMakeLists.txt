cmake_minimum_required(VERSION 3.16)
project(CALOg_NeutronStudy)

# ── Find Geant4 (with UI and VIS optional) ──────────────────────────────────
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

# ── Sources ──────────────────────────────────────────────────────────────────
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(calog_sim calog_sim.cc ${SOURCES})
target_link_libraries(calog_sim ${Geant4_LIBRARIES})

# ── Copy macros to build dir ─────────────────────────────────────────────────
configure_file(run.mac  ${CMAKE_CURRENT_BINARY_DIR}/run.mac  COPYONLY)

