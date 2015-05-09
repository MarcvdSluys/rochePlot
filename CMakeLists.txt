##
##  CMakeFile for rochePlot
##  MvdS, 2015-05-06
##
##  Copyright 1993-2015 Frank Verbunt, Marc van der Sluys
##  
##  This file is part of the rochePlot package,
##  see: http://rocheplot.sf.net/
##  
##  RochePlot is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##  
##  RochePlot is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##  
##  You should have received a copy of the GNU General Public License
##  along with rochePlot.  If not, see <http://www.gnu.org/licenses/>.
##
##  To install this software on your system, see the INSTALL file
##  in the root or doc/ directory of this package.
##


cmake_minimum_required( VERSION 2.8 )

# Set build type. Do this *before* we set the project name:
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo Profile."
    FORCE )
endif( NOT CMAKE_BUILD_TYPE )
set( CMAKE_CONFIGURATION_TYPES "${CMAKE_BUILD_TYPE}" CACHE INTERNAL "internal" )


# Project name and language:
project( rochePlot Fortran )


# Search in the CMake/ directory for CMake modules:
list( APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake )


# Various compile/optimisation options that we may want to enable:
include( SetCompileOptions )


# Get the compiler name (for compiler flags and to find libraries):
get_filename_component( Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME )



# Find dependencies:
find_package( PGPLOT REQUIRED )
find_package( LibSUFR REQUIRED )
set( INCLUDE_FLAGS "-I${LibSUFR_INCLUDES}" )  # will be transferred to CompilerFlags


# Set FORTRAN compiler flags:
include( CompilerFlags_Fortran )


# Create the file code_version.f90, which contains the code version number/hash and date:
if( NOT EXISTS ${CMAKE_SOURCE_DIR}/src/code_version.f90 OR CREATE_VERSION )
  # Code version generator:
  add_custom_command(
    OUTPUT ${CMAKE_SOURCE_DIR}/src/code_version.f90
    COMMAND cd $(CMAKE_SOURCE_DIR)
    COMMAND . ${CMAKE_SOURCE_DIR}/code_version.sh $(CMAKE_SOURCE_DIR) src/code_version.f90 ${Fortran_COMPILER_NAME} ${OPT_FLAGS}
    )
  
  # Tell CMake the source won't be available until build time:
  set_source_files_properties( ${CMAKE_SOURCE_DIR}/src/code_version.f90 PROPERTIES GENERATED 1 )
endif( NOT EXISTS ${CMAKE_SOURCE_DIR}/src/code_version.f90 OR CREATE_VERSION )


# Compile and link the code:
add_executable( rochePlot  src/rochePlot.f90 src/code_version.f90 )
target_link_libraries( rochePlot ${PGPLOT_LIBRARIES} ${LibSUFR_LIBRARIES} )

# Install the binary and example input file:
install( TARGETS  rochePlot   RUNTIME DESTINATION bin )
install( FILES    share/rochePlot.dat DESTINATION share/rochePlot )
