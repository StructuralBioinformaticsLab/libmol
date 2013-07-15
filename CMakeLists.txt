if(NOT EMBED_LIBMOL)
  cmake_minimum_required(VERSION 2.8)
endif()

project(libmol)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules/")

set(libmol_version_major 0)
set(libmol_version_minor 0)
set(libmol_version_patch 6)
set(libmol_version
  "${libmol_version_major}.${libmol_version_minor}.${libmol_version_patch}")

set(SOURCES
  mol.0.0.6/atom.c
  mol.0.0.6/atom_group.c
  mol.0.0.6/_atom_group_copy_from_deprecated.c
  mol.0.0.6/benergy.c
  mol.0.0.6/bond.c
  mol.0.0.6/compare.c
  mol.0.0.6/energy.c
  mol.0.0.6/gbsa.c
  mol.0.0.6/hbond.c
  mol.0.0.6/icharmm.c
  mol.0.0.6/init.c
  mol.0.0.6/io.c
  mol.0.0.6/json.c
  mol.0.0.6/lbfgs.c
  mol.0.0.6/matrix.c
  mol.0.0.6/mem.c
  mol.0.0.6/minimize.c
  mol.0.0.6/mol2.c
  mol.0.0.6/ms.c
  mol.0.0.6/myhelpers.c
  mol.0.0.6/nbenergy.c
  mol.0.0.6/octree.c
  mol.0.0.6/pdb.c
  mol.0.0.6/potential.c
  mol.0.0.6/prms.c
  mol.0.0.6/protein.c
  mol.0.0.6/rmsd.c
  mol.0.0.6/rotamer.c
  mol.0.0.6/sasa.c
  mol.0.0.6/sdf.c
  mol.0.0.6/shield.c
  mol.0.0.6/subag.c
  mol.0.0.6/version.c)

add_library(mol.${libmol_version} ${SOURCES})

include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/mol.0.0.6/")

find_package(Git)
if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --always HEAD
    COMMAND sed -e "s/-/./g"
    OUTPUT_VARIABLE GIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
  set(GIT_VERSION "git")
endif()

add_definitions(
  -D _MOL_VERSION_="\\"${libmol_version}\\""
  -D _MOL_INCLUDE_="\\"../mol.${libmol_version}.h\\""
  -D _GIT_VERSION_="\\"${GIT_VERSION}\\""
  -ffast-math
  -Wall
  -W
  -Wshadow
  -Wpointer-arith
  -Wcast-qual
  -std=c99
  -Winline)

if(NOT CMAKE_SYSTEM_PROCESSOR STREQUAL "ppc64")
  set(LIB_INSTALL_DIR "$ENV{HOME}/lib")
  set(HEADER_INSTALL_DIR "$ENV{HOME}/include")
else()
  # assume BG/L
  add_definitions(
    -D _BGL_
    -D _MPI_)
  include_directories(
    "/bgl/BlueLIght/ppcfloor/bglsys/include"
    "$ENV{HOME}/usr/include")
  link_directories(
    "/bgl/BlueLight/ppcfloor/bglsys/lib"
    "$ENV{HOME}/usr/lib")

  set(LIB_INSTALL_DIR "$ENV{HOME}/usr/lib")
  set(HEADER_INSTALL_DIR "$ENV{HOME}/usr/include")
endif()

if(APPLE)
  add_definitions(
    -D _DARWIN_)
  if(_CURRENT_OSX_VERSION VERSION_EQUAL 10.6)
    add_definitions(
      -D _DARWIN_SNOW_LEOPARD_)
  endif()
endif()

# Testing
find_package(Check)

if(CHECK_FOUND)
  enable_testing()
  add_subdirectory(tests)
else()
  message(Could not find Check for testing)
endif()
