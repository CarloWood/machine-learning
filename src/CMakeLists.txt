cmake_minimum_required(VERSION 3.12...3.27)

find_package(Threads REQUIRED)
find_package(TensorflowCC REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(OPENMP REQUIRED openmp)
pkg_check_modules(CAIRO REQUIRED cairo)
pkg_check_modules(X11 REQUIRED x11)

include(AICxxProject)

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
target_link_libraries(logistic gnuplot_i::gnuplot_i)

add_executable(likeliness likeliness.cxx)
target_link_libraries(likeliness gnuplot_i::gnuplot_i)

add_executable(residual_loss residual_loss.cxx)
target_link_libraries(residual_loss gnuplot_i::gnuplot_i ${AICXX_OBJECTS_LIST} ${CAIRO_LDFLAGS} ${X11_LDFLAGS})

add_executable(cairo_window cairo_window.cxx)
target_compile_options(cairo_window PRIVATE ${CAIRO_CFLAGS} ${X11_CFLAGS})
target_link_libraries(cairo_window ${AICXX_OBJECTS_LIST} ${CAIRO_LDFLAGS} ${X11_LDFLAGS})

add_executable(quadratic_bezier quadratic_bezier.cxx)
target_compile_options(quadratic_bezier PRIVATE ${CAIRO_CFLAGS} ${X11_CFLAGS})
target_link_libraries(quadratic_bezier ${AICXX_OBJECTS_LIST} ${CAIRO_LDFLAGS} ${X11_LDFLAGS})

add_subdirectory(tfdebug)
add_subdirectory(intro)
add_subdirectory(ops)
add_subdirectory(tests)
