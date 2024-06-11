#! /bin/sh

# Run the initial command
echo -n "Running as: "; whoami
echo "Current Working directory: \"$PWD\""

echo "DISPLAY = \"${DISPLAY}\""
echo "GITACHE_ROOT = \"${GITACHE_ROOT}\""
echo "AUTOGEN_CMAKE_ONLY = ${AUTOGEN_CMAKE_ONLY}"
echo "TensorflowCC_DIR = \"${TensorflowCC_DIR}\""
echo "TENSORFLOW_PREFIX = \"${TENSORFLOW_PREFIX}\""

if [ -z "$DISPLAY" ]; then
  echo "Warning: DISPLAY is not set!"
fi
if [ -z "$GITACHE_ROOT" ]; then
  export GITACHE_ROOT="$HOME/data/gitache"
fi
if [ -z "$AUTOGEN_CMAKE_ONLY" ]; then
  export AUTOGEN_CMAKE_ONLY=1
fi
if [ -z "$TensorflowCC_DIR" ]; then
  export TensorflowCC_DIR="$HOME/data/machine-learning/cmake"
fi
if [ -z "$TENSORFLOW_PREFIX" ]; then
  export TENSORFLOW_PREFIX=/usr
fi

if [ ! -d "$GITACHE_ROOT" ]; then
  mkdir $GITACHE_ROOT
fi

if ! command -v git &> /dev/null; then
  sudo pacman --noconfirm -Sy git clang make pkgconfig cmake boost pango cairo eigen tensorflow-opt neovim ttf-dejavu

  echo "Now enter the container with the command:"
  echo "docker exec -it machine-learning su --login --whitelist-environment="DISPLAY,GITACHE_ROOT,AUTOGEN_CMAKE_ONLY,TensorflowCC_DIR,TENSORFLOW_PREFIX" archuser"
  echo "Then download and configure the project:"
  echo "$ cd data"
  echo "$ git clone --recursive https://github.com/CarloWood/machine-learning.git"
  echo "$ cd machine-learning"
  echo "$ mkdir build"
  echo "$ cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DCMAKE_VERBOSE_MAKEFILE=ON -DEnableCairoWindowTests:BOOL=ON"
  echo "$ To build the project run:"
  echo "$ cmake --build build --config Debug --parallel $(nproc --all)"
  echo "$ cd build/cairowindow/tests"
  echo "$ ./quadratic_bezier"

  # Hang here to preserve everything above. Use the above `docker exec -it machine-learning ...` command to access the container.
  tail -f /dev/null
fi

if [ ! -d machine-learning ]; then
  git clone --recursive https://github.com/CarloWood/machine-learning.git
fi

cd machine-learning

if [ ! -d build ]; then
  mkdir build
fi

cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DCMAKE_VERBOSE_MAKEFILE=ON -DEnableCairoWindowTests:BOOL=ON
cmake --build build --config Debug --parallel $(nproc --all)

if [ $? -eq 0 ]; then
  tail -f /dev/null
else
  echo "An error occurred: leaving docker container."
fi
