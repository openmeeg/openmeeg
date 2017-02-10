if [[ "$USE_PROJECT" == "0" ]]; then
  cd OpenMEEG
fi

mkdir build
cd build

# XXX : BUILD_SHARED should be used to set global defaults

function install_matio {  # Install MATIO
  wget https://github.com/openmeeg/matio-openmeeg/archive/master.zip
  unzip master.zip
  cd matio-openmeeg-master
  mkdir build
  cd build
  cmake -DMAT73:BOOL=ON ..
  make
  sudo make install
  cd ../../
}

if [[ $USE_PYTHON == "1" ]]; then
  ENABLE_PYTHON=ON
else
  ENABLE_PYTHON=OFF
fi

if [[ $USE_VTK == "1" ]]; then
  USE_VTK=ON
else
  USE_VTK=OFF
fi

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
  if [[ "$USE_PROJECT" == "0" ]]; then
      # Build OpenMEEG and use openblas from homebrew
      # or use Atlas which maps to vecLib or Accelerate frameworks
      cmake \
      -DBUILD_SHARED:BOOL=ON \
      -DBUILD_DOCUMENTATION:BOOL=OFF \
      -DBUILD_TESTING:BOOL=ON \
      -DENABLE_PYTHON:BOOL=${ENABLE_PYTHON} \
      -DENABLE_PACKAGING:BOOL=ON \
      -DUSE_VTK:BOOL=${USE_VTK} \
      -DBLASLAPACK_IMPLEMENTATION:BOOL=${BLASLAPACK_IMPLEMENTATION} \
      -DCMAKE_SKIP_RPATH:BOOL=OFF \
      ..
  else
      cmake \
      -DBUILD_DOCUMENTATION:BOOL=OFF \
      -DBUILD_TESTING:BOOL=ON \
      -DENABLE_PYTHON:BOOL=${ENABLE_PYTHON} \
      -DENABLE_PACKAGING:BOOL=ON \
      -DUSE_VTK:BOOL=${USE_VTK} \
      -DBLASLAPACK_IMPLEMENTATION:BOOL=${BLASLAPACK_IMPLEMENTATION} \
      -DUSE_SYSTEM_matio:BOOL=OFF \
      -DUSE_SYSTEM_hdf5:BOOL=OFF \
      -DCMAKE_SKIP_RPATH:BOOL=OFF \
      ..
  fi
else
  if [[ "$USE_SYSTEM" == "1" ]]; then
    if [[ "$USE_PROJECT" == "0" ]]; then
      install_matio

        # Build OpenMEEG with ATLAS or openblas
        cmake \
          -DBUILD_SHARED:BOOL=ON \
          -DBUILD_DOCUMENTATION:BOOL=OFF \
          -DBUILD_TESTING:BOOL=ON \
          -DENABLE_PYTHON:BOOL=${ENABLE_PYTHON} \
          -DENABLE_PACKAGING:BOOL=ON \
          -DUSE_VTK:BOOL=${USE_VTK} \
          -DBLASLAPACK_IMPLEMENTATION:BOOL=${BLASLAPACK_IMPLEMENTATION} \
          -DCMAKE_SKIP_RPATH:BOOL=OFF \
          ..
    else
      cmake \
          -DBUILD_DOCUMENTATION:BOOL=OFF \
          -DBUILD_TESTING:BOOL=ON \
          -DENABLE_PYTHON:BOOL=${ENABLE_PYTHON} \
          -DENABLE_PACKAGING:BOOL=ON \
          -DUSE_VTK:BOOL=${USE_VTK} \
          -DBLASLAPACK_IMPLEMENTATION:BOOL=${BLASLAPACK_IMPLEMENTATION} \
          -DUSE_SYSTEM_matio:BOOL=OFF \
          -DUSE_SYSTEM_hdf5:BOOL=ON \
          -DUSE_SYSTEM_zlib:BOOL=ON \
          -DCMAKE_SKIP_RPATH:BOOL=OFF \
          ..
      fi
  else
    cmake \
          -DBUILD_DOCUMENTATION:BOOL=OFF \
          -DBUILD_TESTING:BOOL=ON \
          -DENABLE_PYTHON:BOOL=${ENABLE_PYTHON} \
          -DENABLE_PACKAGING:BOOL=ON \
          -DBLASLAPACK_IMPLEMENTATION:BOOL=${BLASLAPACK_IMPLEMENTATION} \
          -DUSE_VTK:BOOL=${USE_VTK} \
          -DUSE_SYSTEM_matio:BOOL=OFF \
          -DUSE_SYSTEM_hdf5:BOOL=OFF \
          -DUSE_SYSTEM_zlib:BOOL=OFF \
          -DCMAKE_SKIP_RPATH:BOOL=OFF \
          ..
  fi
fi

make

ctest -V
