#==========================================================================
#  AIDA Detector description implementation 
#--------------------------------------------------------------------------
# Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
# All rights reserved.
#
# For the licensing terms see $DD4hepINSTALL/LICENSE.
# For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
#
#==========================================================================
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project("JunoDetector")
include_directories(include)

IF(NOT TARGET DD4hep::DDCore)
  find_package ( DD4hep REQUIRED )
  include ( ${DD4hep_DIR}/cmake/DD4hep.cmake )
  include ( ${DD4hep_DIR}/cmake/DD4hepBuild.cmake )
  dd4hep_configure_output()
ENDIF()

dd4hep_set_compiler_flags()
#==========================================================================
if(TARGET XercesC::XercesC)
  SET(OPT_XERCESC XercesC::XercesC)
endif()

dd4hep_add_plugin(mmDDDetectors SOURCES src/*.cpp
  USES DD4hep::DDCore DD4hep::DDRec ROOT::Core ROOT::Geom ROOT::GenVector ${OPT_XERCESC}
  )
install(TARGETS mmDDDetectors LIBRARY DESTINATION lib)
install(DIRECTORY compact DESTINATION examples/ClientTests)
#-----------------------------------------------------------------------------------
dd4hep_configure_scripts(mmDDDetectors DEFAULT_SETUP WITH_TESTS)

#***  Testing  *********************************************************************
foreach (test simpletttest_ILD)
##  dd4hep_add_test_reg( DDDetectors_converter_gdml
##    COMMAND    "${CMAKE_INSTALL_PREFIX}/bin/run_test_DDDetectors.sh"
##    EXEC_ARGS  geoConverter -compact2gdml -input file:${CMAKE_CURRENT_SOURCE_DIR}/compact/${test}.xml -output file:${test}.gdml
##    REGEX_PASS " Successfully extracted GDML to")
  foreach( type description vis )
    dd4hep_add_test_reg( DDDetectors_converter_${type}
      COMMAND    "${CMAKE_INSTALL_PREFIX}/bin/run_test_DDDetectors.sh"
      EXEC_ARGS  geoConverter -compact2${type} -input file:${CMAKE_CURRENT_SOURCE_DIR}/compact/${test}.xml -output file:${test}.${type}
      REGEX_PASS " Handled [1-9][0-9]* volumes")
  endforeach(type)
endforeach(test)


add_subdirectory(dd4hep2FBX)
# set(exename dd4hep2FBX)
# include_directories(${PROJECT_SOURCE_DIR}/dd4hep2FBX)
# add_executable(${exename} test.cc)
# target_link_libraries(${exename} DD4hep::DDCore DD4hep::DDCond 
#         DD4hep::DDAlign ROOT::Core ROOT::Geom ROOT::GenVector 
#         libdd4hep2FBX)


#==========================================================================
#
#  AIDA Detector description implementation 
#--------------------------------------------------------------------------
# Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
# All rights reserved.
#
# For the licensing terms see $DD4hepINSTALL/LICENSE.
# For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
#
#==========================================================================
#--------------------------------------------------------------------------
#  DDDetectors
#
#  Author:  M.Frank
#  Date:    13/11/2014
#
#--------------------------------------------------------------------------
# SET_PROPERTY(DIRECTORY . PROPERTY PACKAGE_NAME DDDetectors)

# dd4hep_add_plugin(DDDetectors 
#   SOURCES src/*.cpp
#   USES DD4hep::DDRec
#   INCLUDES $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
# )

# set_target_properties(DDDetectors PROPERTIES VERSION ${DD4hep_VERSION} SOVERSION ${DD4hep_SOVERSION})
# install(TARGETS DDDetectors EXPORT DD4hep
#   LIBRARY DESTINATION lib
#   ARCHIVE DESTINATION lib
#   )
# install(DIRECTORY compact DESTINATION DDDetectors)

