cmake_minimum_required(VERSION 3.12...3.27)

# Test written by ChatGPT.
add_test(set_inter_op_parallelism_threads set_inter_op_parallelism_threads)

# Test for calculating intersection point of two lines (Line::intersection_with).
add_executable(intersecting_lines intersecting_lines.cxx)
target_compile_options(intersecting_lines PRIVATE ${CAIRO_CFLAGS} ${X11_CFLAGS})
target_link_libraries(intersecting_lines ${AICXX_OBJECTS_LIST} ${CAIRO_LDFLAGS} ${X11_LDFLAGS})
