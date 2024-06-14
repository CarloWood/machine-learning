This is an ongoing effort to do Machine Learning in pure C++,
leveraging Tensorflow to make it run on the GPU.

To compile this project you need linux with TensorflowCC;
since the latter is rather hard to compile yourself, I recommend
to use Arch (or another distribution if that provides a pre-packed
tensorflow_cc).

USING DOCKER
============

As an alternative you can work on this project in a docker
container; which even works on Windows! Start the docker
container from wsl in that case (installing the project
directly in wsl doesn't work very well it seems).

Read https://github.com/CarloWood/machine-learning/tree/master/docker#readme
for an explanation of how to construct the container and run it.

Once inside the docker contain (you are then in the directory
/home/archuser/machine-learning), run:

```
                   docker:~/machine-learning> git clone --recursive https://github.com/CarloWood/machine-learning.git
                   docker:~/machine-learning> cd machine-learning
  docker:~/machine-learning/machine-learning> cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DCMAKE_VERBOSE_MAKEFILE=ON -DEnableCairoWindowTests:BOOL=ON
  docker:~/machine-learning/machine-learning> cmake --build build --config Debug --parallel $(nproc --all)
  docker:~/machine-learning/machine-learning> cd build/cairowindow/tests
  docker:~/machine-learning/machine-learning> ./quadratic_bezier
```

COMPILING ON A LINUX HOST
=========================

This works the same as in docker, except that in this case you have to set
up the build environment yourself.

On Archlinux you need to following packages:

  git clang make pkgconfig cmake boost pango cairo eigen tensorflow-opt ttf-dejavu ccache

The last two are not absolutely required, but you really want to install them.


