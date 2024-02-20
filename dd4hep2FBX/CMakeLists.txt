set(exename dd4hep2FBX)
add_executable(${exename} dd4hep2FBXWriter.cc)

target_link_libraries(${exename} DD4hep::DDCore DD4hep::DDCond DD4hep::DDAlign ROOT::Core ROOT::Geom ROOT::GenVector )

install(TARGETS ${exename}
  # IMPORTANT: Add the bar executable to the "export-set"
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)