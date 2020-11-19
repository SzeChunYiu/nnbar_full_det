cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(nnbar-calo-sim)

option(WITH_GEANT4_UIVIS "Build nnbar-calo-sim with Geant4 UI and Vis drivers" ON)

if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#

include_directories(${PROJECT_SOURCE_DIR}/include )

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*)

add_executable(nnbar-calo-sim nnbar-calo-sim.cc ${sources} ${headers})
target_link_libraries(nnbar-calo-sim ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

set(SCRIPTS
  nnbar-calo-sim.in
  gui.mac
  init_vis.mac
  vis.mac
  setup.file
  )

foreach(_script ${SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS nnbar-calo-sim DESTINATION bin)
