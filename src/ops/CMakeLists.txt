cmake_minimum_required(VERSION 3.12...3.27)

add_test(switchmerge SwitchMerge)
target_include_directories(switchmerge PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(switchmerge tfdebug::tfdebug)
