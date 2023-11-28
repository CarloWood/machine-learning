cmake_minimum_required(VERSION 3.12...3.27)

# Here simply call the macro for every sample file (the one that contains the main).
# These tests comes from https://github.com/ksachdeva/tensorflow-cc-examples/tree/master/examples/1_Introduction/src/.
add_test(intro-hello-world hello-world)
add_test(intro-basic-operations basic-operations)
add_test(intro-naming-ops naming-ops)
add_test(intro-basic-tensors basic-tensors)
add_test(intro-basic-graphdef basic-graphdef)
add_test(intro-basic-io basic-io)
