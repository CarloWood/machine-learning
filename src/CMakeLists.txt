cmake_minimum_required(VERSION 3.12...3.27)

find_package(Threads REQUIRED)
find_package(TensorflowCC REQUIRED)
#find_package(CUDA)

find_package(PkgConfig REQUIRED)
pkg_check_modules(OPENMP REQUIRED openmp)

# Macro tp help create the executables.
macro(add_test app_name main_file)
  add_executable(${app_name} ${main_file}.cxx)
  target_compile_options(${app_name} PRIVATE ${OPENMP_CFLAGS})
  target_link_libraries(${app_name} TensorflowCC::TensorflowCC ${OPENMP_LDFLAGS} Threads::Threads)
endmacro()

# The example used by https://github.com/FloopCZ/tensorflow_cc/tree/master/example/.
add_test(floopcz_example floopcz_example)

# Test written by ChatGPT.
add_test(set_inter_op_parallelism_threads set_inter_op_parallelism_threads)

# Here simply call the macro for every sample file (the one that contains the main).
# These tests comes from https://github.com/ksachdeva/tensorflow-cc-examples/tree/master/examples/1_Introduction/src/.
add_test(intro-hello-world hello-world)
add_test(intro-basic-operations basic-operations)
add_test(intro-naming-ops naming-ops)
add_test(intro-basic-tensors basic-tensors)
add_test(intro-basic-graphdef basic-graphdef)
add_test(intro-basic-io basic-io)
