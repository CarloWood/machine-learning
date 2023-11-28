cmake_minimum_required(VERSION 3.12...3.27)

project(tfdebug
  DESCRIPTION "Debugging support for tensorflow programs."
  LANGUAGES CXX
)

find_package(TensorflowCC REQUIRED)

# This project is an OBJECT-library, used by other git submodules and the main project.
add_library(tfdebug_ObjLib OBJECT)
target_link_libraries(tfdebug_ObjLib TensorflowCC::TensorflowCC)

# The list of source files.
target_sources(tfdebug_ObjLib
  PRIVATE
    "dump_graph.cxx"
    "DebugLayer.h"
    "DebugShow.h"
    "DebugTensor.h"
)

target_include_directories(tfdebug_ObjLib
  INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

# Create an ALIAS target.
add_library(tfdebug::tfdebug ALIAS tfdebug_ObjLib)
