cmake_minimum_required(VERSION 3.24)

project(
  libfci
  VERSION 0.1.0
  DESCRIPTION "Compact determinant full-CI solver library"
  HOMEPAGE_URL "https://github.com/DCM-Uni-Paderborn/libfci"
  LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(LIBFCI_USE_OPENMP "Enable OpenMP parallel regions when available" ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "RelWithDebInfo"
      CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release
                                               RelWithDebInfo)
endif()

find_package(LAPACK REQUIRED)

add_library(
  fci
  src/BinomialCoeff.cpp
  src/CxAlgebra.cpp
  src/CxAssertFail.cpp
  src/CxMemoryStack.cpp
  src/CxStorageDevice.cpp
  src/Fci.cpp
  src/FciCsf.cpp
  src/FciDump.cpp
  src/FciOps.cpp
  src/FciSubspace.cpp
  src/FciVec.cpp)

add_library(libfci::fci ALIAS fci)

target_compile_features(fci PRIVATE cxx_std_11)
target_include_directories(
  fci PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
             $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(fci PRIVATE LAPACK::LAPACK)

set_target_properties(
  fci
  PROPERTIES VERSION ${PROJECT_VERSION}
             SOVERSION ${PROJECT_VERSION_MAJOR}
             POSITION_INDEPENDENT_CODE ON)

if(LIBFCI_USE_OPENMP)
  find_package(OpenMP QUIET COMPONENTS CXX)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(fci PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()

set(LIBFCI_PC_PRIVATE_LIBS "${LAPACK_LIBRARIES}")
string(REPLACE ";" " " LIBFCI_PC_PRIVATE_LIBS "${LIBFCI_PC_PRIVATE_LIBS}")

install(TARGETS fci EXPORT libfci-targets)
install(FILES include/libfci.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(
  EXPORT libfci-targets
  NAMESPACE libfci::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libfci)

configure_package_config_file(
  cmake/libfci-config.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/libfci-config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libfci)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/libfci-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libfci-config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/libfci-config-version.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libfci)

configure_file(cmake/libfci.pc.in ${CMAKE_CURRENT_BINARY_DIR}/libfci.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libfci.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
