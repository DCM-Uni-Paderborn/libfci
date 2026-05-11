cmake_minimum_required(VERSION 3.16)

project(
  libfci
  VERSION 0.1.0
  DESCRIPTION "Compact determinant full-CI solver library"
  HOMEPAGE_URL "https://github.com/DCM-Uni-Paderborn/libfci"
  LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(LIBFCI_USE_OPENMP "Enable OpenMP parallel regions when available" ON)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(LIBFCI_SOURCES
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

add_library(fci ${LIBFCI_SOURCES})
add_library(libfci::fci ALIAS fci)

target_compile_features(fci PUBLIC cxx_std_11)
target_include_directories(
  fci
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(fci PUBLIC LAPACK::LAPACK BLAS::BLAS)

set_target_properties(
  fci
  PROPERTIES OUTPUT_NAME fci
             VERSION ${PROJECT_VERSION}
             SOVERSION ${PROJECT_VERSION_MAJOR}
             POSITION_INDEPENDENT_CODE ON)

if(LIBFCI_USE_OPENMP)
  find_package(OpenMP QUIET COMPONENTS CXX)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(fci PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()

set(LIBFCI_PC_PRIVATE_LIBS "${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}")
string(REPLACE ";" " " LIBFCI_PC_PRIVATE_LIBS "${LIBFCI_PC_PRIVATE_LIBS}")

install(
  TARGETS fci
  EXPORT libfciTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES include/libfci.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(
  EXPORT libfciTargets
  NAMESPACE libfci::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libfci)

configure_package_config_file(
  cmake/libfciConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/libfciConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libfci)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/libfciConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libfciConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/libfciConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libfci)

configure_file(libfci.pc.in libfci.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libfci.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
