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
    # sudo apt-get install doxygen libatlas-dev libatlas-base-dev libblas-dev liblapack-dev # python-numpy swig python-dev libvtk5-dev libtiff4-dev
    # sudo apt-get install doxygen libblas-dev liblapack-dev # python-numpy swig python-dev libvtk5-dev libtiff4-dev
    # sudo apt-get install doxygen python-numpy swig python-dev
    sudo apt-get install doxygen sshpass
    wget https://s3.amazonaws.com/biibinaries/thirdparty/cmake-3.0.2-Linux-64.tar.gz
    tar -xzf cmake-3.0.2-Linux-64.tar.gz
    sudo cp -fR cmake-3.0.2-Linux-64/* /usr
    rm -rf cmake-3.0.2-Linux-64
    rm cmake-3.0.2-Linux-64.tar.gz
fi

mkdir build
cd build

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
  cmake \
      -DATLAS_INCLUDE_PATH:PATH=/usr/include/atlas \
      -DBUILD_SHARED:BOOL=ON \
      -DBUILD_DOCUMENTATION:BOOL=OFF \
      -DBUILD_TESTING:BOOL=ON \
      -DENABLE_PYTHON:BOOL=OFF \
      -DENABLE_PACKAGING:BOOL=ON \
      -DUSE_VTK:BOOL=OFF \
      -DUSE_ATLAS:BOOL=ON \
      -DUSE_SYSTEM_MATIO:BOOL=OFF \
      -DUSE_SYSTEM_hdf5:BOOL=OFF \
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
      -DUSE_ATLAS:BOOL=OFF \
      -DUSE_SYSTEM_MATIO:BOOL=OFF \
      -DUSE_SYSTEM_hdf5:BOOL=OFF \
      -DCMAKE_SKIP_RPATH:BOOL=OFF \
      ..
fi

make
