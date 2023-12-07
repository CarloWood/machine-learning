cmake_minimum_required(VERSION 3.12...3.27)

find_package(Threads REQUIRED)
find_package(TensorflowCC REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(OPENMP REQUIRED openmp)

# Macro to help create the executables.
macro(add_test app_name main_file)
  add_executable(${app_name} ${main_file}.cxx)
  target_compile_options(${app_name} PRIVATE ${OPENMP_CFLAGS})
  target_link_libraries(${app_name} TensorflowCC::TensorflowCC ${OPENMP_LDFLAGS} Threads::Threads)
endmacro()

# My stuff.
add_test(neural_network neural_network)
target_link_libraries(neural_network tfdebug::tfdebug gnuplot_i::gnuplot_i)

add_executable(logistic logistic.cxx)
target_include_directories(logistic PRIVATE "${CMAKE_SOURCE_DIR}/stats/include")
target_link_libraries(logistic gnuplot_i::gnuplot_i)

add_executable(likeliness likeliness.cxx)
target_include_directories(likeliness PRIVATE "${CMAKE_SOURCE_DIR}/stats/include")
target_link_libraries(likeliness gnuplot_i::gnuplot_i)

add_subdirectory(tfdebug)
add_subdirectory(intro)
add_subdirectory(ops)
add_subdirectory(tests)
