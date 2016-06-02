if [[ $TRAVIS_OS_NAME == 'osx' ]]; then

    # Install some custom requirements on OS X
    $CXX --version

else
    # # g++4.8.1
    if [ "$CXX" = "g++" ]; then sudo apt-get install -qq g++-4.8; fi
    if [ "$CXX" = "g++" ]; then export CXX="g++-4.8"; fi

    # clang 3.4
    if [ "$CXX" == "clang++" ]; then sudo apt-get install --allow-unauthenticated -qq clang-3.4; fi
    if [ "$CXX" == "clang++" ]; then export CXX="clang++-3.4"; fi

    sudo apt-get update
    if [[ $USE_SYSTEM == "1" ]]; then
      sudo apt-get install libatlas-dev libatlas-base-dev libblas-dev liblapack-dev libhdf5-dev
      # python-numpy swig python-dev libvtk5-dev libtiff4-dev doxygen
    fi
    # sudo apt-get install libblas-dev libatlas-base-dev liblapack-dev
    # sudo apt-get install libvtk5-dev libtiff4-dev
    # sudo apt-get install python-numpy swig python-dev
    # sudo apt-get install doxygen
    wget https://s3.amazonaws.com/biibinaries/thirdparty/cmake-3.0.2-Linux-64.tar.gz
    tar -xzf cmake-3.0.2-Linux-64.tar.gz
    sudo cp -fR cmake-3.0.2-Linux-64/* /usr
    rm -rf cmake-3.0.2-Linux-64
    rm cmake-3.0.2-Linux-64.tar.gz
fi

if [[ "$USE_PROJECT" == "0" ]]; then
  cd OpenMEEG
fi

mkdir build
cd build

# XXX : BUILD_SHARED should be used to set global defaults

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
  if [[ "$USE_PROJECT" == "0" ]]; then
      # Install MATIO
      wget https://github.com/openmeeg/matio-openmeeg/archive/master.zip
      unzip master.zip
      cd matio-openmeeg-master
      mkdir build
      cd build
      cmake -DMAT73:BOOL=$MAT73 ..
      make
      make install
      cd ../../

      # Build OpenMEEG
      cmake \
      -DBUILD_SHARED:BOOL=ON \
      -DBUILD_DOCUMENTATION:BOOL=OFF \
      -DBUILD_TESTING:BOOL=ON \
      -DENABLE_PYTHON:BOOL=OFF \
      -DENABLE_PACKAGING:BOOL=ON \
      -DUSE_VTK:BOOL=OFF \
      -DUSE_ATLAS:BOOL=ON \
      -DCMAKE_SKIP_RPATH:BOOL=OFF \
      ..
  else
      cmake \
      -DATLAS_INCLUDE_PATH:PATH=/usr/include/atlas \
      -DBUILD_SHARED:BOOL=ON \
      -DBUILD_DOCUMENTATION:BOOL=OFF \
      -DBUILD_TESTING:BOOL=ON \
      -DENABLE_PYTHON:BOOL=OFF \
      -DENABLE_PACKAGING:BOOL=ON \
      -DUSE_VTK:BOOL=OFF \
      -DUSE_ATLAS:BOOL=ON \
      -DUSE_SYSTEM_matio:BOOL=OFF \
      -DUSE_SYSTEM_hdf5:BOOL=OFF \
      -DCMAKE_SKIP_RPATH:BOOL=OFF \
      ..
  fi
else
  if [[ "$USE_SYSTEM" == "1" ]]; then
    cmake \
        -DATLAS_INCLUDE_PATH:PATH=/usr/include/atlas \
        -DBUILD_DOCUMENTATION:BOOL=OFF \
        -DBUILD_TESTING:BOOL=ON \
        -DENABLE_PYTHON:BOOL=OFF \
        -DENABLE_PACKAGING:BOOL=ON \
        -DUSE_VTK:BOOL=OFF \
        -DUSE_ATLAS:BOOL=OFF \
        -DUSE_SYSTEM_matio:BOOL=ON \
        -DUSE_SYSTEM_hdf5:BOOL=ON \
        -DUSE_SYSTEM_zlib:BOOL=ON \
        -DCMAKE_SKIP_RPATH:BOOL=OFF \
        ..
  else
    cmake \
          -DATLAS_INCLUDE_PATH:PATH=/usr/include/atlas \
          -DBUILD_DOCUMENTATION:BOOL=OFF \
          -DBUILD_TESTING:BOOL=ON \
          -DENABLE_PYTHON:BOOL=OFF \
          -DENABLE_PACKAGING:BOOL=ON \
          -DUSE_VTK:BOOL=OFF \
          -DUSE_ATLAS:BOOL=OFF \
          -DUSE_SYSTEM_matio:BOOL=OFF \
          -DUSE_SYSTEM_hdf5:BOOL=OFF \
          -DUSE_SYSTEM_zlib:BOOL=OFF \
          -DCMAKE_SKIP_RPATH:BOOL=OFF \
          ..
  fi
fi

# -DBUILD_SHARED:BOOL=ON \
# -DBUILD_SHARED_LIBS_OpenMEEG:BOOL=ON
# -DBUILD_SHARED_LIBS_hdf5:BOOL=ON -DBUILD_SHARED_LIBS_matio:BOOL=ON
# -DBUILD_SHARED_LIBS_zlib:BOOL=ON

make

ctest -V