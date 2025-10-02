set(CMAKE_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake")

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/MarlinKinfitConfigVersion.cmake
  VERSION ${MarlinKinfit_VERSION}
  COMPATIBILITY SameMajorVersion
)

export(EXPORT MarlinKinfitTargets NAMESPACE MarlinKinfit:: FILE ${PROJECT_BINARY_DIR}/MarlinKinfitTargets.cmake)

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/MarlinKinfitConfig.cmake.in
  ${PROJECT_BINARY_DIR}/MarlinKinfitConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_CMAKEDIR}/${PROJECT_NAME}
)

install(
  FILES
    ${PROJECT_BINARY_DIR}/MarlinKinfitConfig.cmake
    ${PROJECT_BINARY_DIR}/MarlinKinfitConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR}/${PROJECT_NAME}
)
install(EXPORT MarlinKinfitTargets
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR}/${PROJECT_NAME}
  NAMESPACE MarlinKinfit::
)
