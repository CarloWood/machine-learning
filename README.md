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
  docker:~/machine-learning/machine-learning> ./autogen.sh
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

```
  sudo pacman -S git clang make pkgconfig cmake boost pango cairo eigen tensorflow-opt ttf-dejavu ccache
```

The last two are not absolutely required, but you really want to install them.

If you have an NVidia GPU you could use `tensorflow-opt-cuda` instead (but cuda is VERY large
and annoying to download / upgrade; I recommend against installing it if you do not already
have it anyway).

You also need packages that are on github, but not provided by arch (or typically any distribution).
The project automatically clones, configures, compiles and installs those packages using `gitache`.
For `gitache` to work you need to specify a directory that it can use for all of the above; this
directory must be writable by you (and preferably empty when we start).

For example, if you use `/opt/gitache`, you will need:

```
  export GITACHE_ROOT=/opt/gitache
  sudo mkdir $GITACHE_ROOT
  sudo chown $(id -u):$(id -g) $GITACHE_ROOT
```

Of course, you can put also in your home directory if you think the partition of your
home directory is large enough.

Usually it is necessary to run `./autogen.sh` after cloning a repository for the first time.
In order for that script to only set things up for cmake, and not automake, you want to
set the following environment variable:

```
  export AUTOGEN_CMAKE_ONLY=1
```
or just set it while running autogen.sh, like so: `AUTOGEN_CMAKE_ONLY=1 ./autogen.sh`.

You also need to set the following environment variables in order to find TensorflowCC:
```
  export TENSORFLOW_PREFIX=/usr
  export TensorflowCC_DIR=/thefullpath/to/machine-learning/cmake
```
These work for Archlinux - if your tensorflow is installed elsewhere you might
have to change them! For example, if your install already provides `TensorflowCCConfig.cmake` et al,
then `TensorflowCC_DIR` should point there.

Finally add the following file as `~/.libcwdrc` in your home directory:
```
silent = on
channels_default = off

channels_on = warning,debug,notice

gdb = /usr/bin/gdb -x $REPOBASE/.gdbinit $GDBCOMMANDFILE
xterm = konsole --name "attach_gdb" --hide-menubar --nofork --workdir "$PWD" --geometry 165x24-0-0 -e %s
```
or, if you already have that file, make sure that the debug channels `warning`, `debug` and `notice` are on.
Note how this uses `$REPOBASE`. Therefore, in order to make `attach_gdb()` work from within the
program you should also add to your environment:
```
export REPOBASE=/full/path/to/machine-learning
```
Using the directory name of where you cloned this project. Namely, the repository contains a `.gdbinit` which then will be loaded.

After all of this is set up, you can do the same as listed above for docker (`git clone`, etc).
