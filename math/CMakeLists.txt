cmake_minimum_required(VERSION 3.12...3.28)

include(AICxxProject)

#==============================================================================
# DEPENDENCIES
#

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

#==============================================================================
# BUILD OBJECT LIBRARY
#

# This project is an OBJECT-library, used by other git submodules and the main project.
add_library(math_ObjLib OBJECT)

file(GLOB HEADER_FILES "*.h")

# The list of source files.
target_sources(math_ObjLib
  PRIVATE
    "Polynomial.cxx"
    "CubicPolynomial.cxx"
    "AnalyzedCubic.cxx"
    "Point.cxx"
    "Direction.cxx"
    "Line.cxx"
    "LinePiece.cxx"

    ${HEADER_FILES}
)

# Required include search-paths.
get_target_property(CWDS_INTERFACE_INCLUDE_DIRECTORIES AICxx::cwds INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(math_ObjLib
  PUBLIC
    "${CWDS_INTERFACE_INCLUDE_DIRECTORIES}" # For sys.h and debug.h.
)

# Require support for C++20.
target_compile_features(math_ObjLib
  PUBLIC cxx_std_20
)

# Set link dependencies.
target_link_libraries(math_ObjLib
  PRIVATE
    AICxx::utils
    Eigen3::Eigen
)

# Create an ALIAS target.
add_library(AICxx::math ALIAS math_ObjLib)
