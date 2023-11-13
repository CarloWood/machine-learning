include(CMakeFindDependencyMacro)
find_dependency(Threads)
include("${CMAKE_CURRENT_LIST_DIR}/TensorflowCCTargets.cmake")

# Print 'Found' message exactly once.
get_property(_is_defined
  GLOBAL
  PROPERTY FIND_PACKAGE_CONFIG_MESSAGE_PRINTED_TensorflowCC
  DEFINED
)
if(NOT _is_defined)
  set(FIND_PACKAGE_MESSAGE_DETAILS_TensorflowCC CACHED INTERNAL "")
  set_property(
    GLOBAL
    PROPERTY FIND_PACKAGE_CONFIG_MESSAGE_PRINTED_TensorflowCC 1
  )
endif()
find_package_message(
  TensorflowCC
  "Found TensorflowCC: ${TensorflowCC_DIR} (found version \"${TensorflowCC_VERSION}\")" "[${CMAKE_CURRENT_LIST_FILE}][${TensorflowCC_DIR}][${TensorflowCC_VERSION}]"
)
