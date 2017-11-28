function setup_conda_wrap {
    if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
        curl https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh -o miniconda.sh -s
    else
        wget -q http://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh
    fi

    chmod +x miniconda.sh
    ./miniconda.sh -b -p ${HOME}/miniconda_wrap
    export PATH=${HOME}/miniconda_wrap/bin:$PATH
    conda update --yes --quiet conda
    conda create -n wrappingenv --yes pip python=$PYTHON_VERSION
    conda create -n wrappingenv --yes pip python=$PYTHON_VERSION
    source activate wrappingenv
    conda install -y --quiet numpy swig
}

# install MKL if necessary.

if [[ "$BLASLAPACK_IMPLEMENTATION" == "MKL" ]]; then
    cmake -P cmake/InstallMKL.cmake
fi

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then

    # Install some custom requirements on OS X
    brew tap homebrew/science # a lot of cool formulae for scientific tools
    # brew update && brew upgrade

    # Install some custom requirements on OS X

    pkgs=""
    links=""
    if [[ "$USE_PROJECT" == "0" || "$USE_SYSTEM" == "1" ]]; then
        pkgs="${pkgs} hdf5 libmatio"
    fi

    if [[ "$BLASLAPACK_IMPLEMENTATION" == "OpenBLAS" || "$BLASLAPACK_IMPLEMENTATION" == "Auto" ]]; then
        pkgs="${pkgs} hdf5 openblas"
        links="${links} openblas --force" # required as link is not automatic
    fi

    if [[ "$USE_VTK" == "1" && "$ENABLE_PACKAGING" != "1" ]]; then
        pkgs="${pkgs} vtk"
    fi

    if [[ "$USE_OMP" == 1 ]]; then
        pkgs="${pkgs} llvm"
        export OMP_NUM_THREADS=4
        export LDFLAGS="-L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib"
        export DYLD_LIBRARY_PATH="/usr/local/opt/llvm/lib:$DYLD_LIBRARY_PATH"
        export CC="/usr/local/opt/llvm/bin/clang"
        export CXX="/usr/local/opt/llvm/bin/clang++"
        export CPPFLAGS="-I/usr/local/opt/llvm/include -fopenmp $CPPFLAGS"
        export CFLAGS="-I/usr/local/opt/llvm/include -fopenmp $CFLAGS"
    fi

    if [[ "$USE_CGAL" == 1 ]]; then
        pkgs="${pkgs} cgal"
    fi

    if [[ "$BUILD_DOCUMENTATION" == "ON" ]]; then
        pkgs="${pkgs} Doxygen Graphviz" # For building documentation
    fi

    if [[ "${pkgs}" != "" ]]; then
        brew install ${pkgs}
    fi

    if [[ "${links}" != "" ]]; then
        brew link ${links}
    fi
else
    # Install some custom requirements on Linux
    export CXX="g++"; 

    # clang
    if [ "$CXX" == "clang++" ]; then
        export CXX="clang++";
    fi

    pkgs=""
    if [[ "$USE_PROJECT" == "0" || "$USE_SYSTEM" == "1" ]]; then
        pkgs="${pkgs} libhdf5-serial-dev libmatio-dev"
    fi

    if [[ "$USE_CGAL" == 1 ]]; then
        pkgs="${pkgs} libcgal-dev"
    fi

    if [[ "$USE_GIFTI" == 1 ]]; then
        pkgs="${pkgs} libnifti-dev libgiftiio-dev"
    fi

    if [[ "$BLASLAPACK_IMPLEMENTATION" == "Atlas" ]]; then
        pkgs="${pkgs} libatlas-dev libatlas-base-dev"
    elif [[ "$BLASLAPACK_IMPLEMENTATION" == "LAPACK" ]]; then
        if [[ "$USE_PROJECT" == "0" || "$USE_SYSTEM" == "1" ]]; then
            pkgs="${pkgs} liblapack-dev libblas-dev"
        fi
    elif [[ "$BLASLAPACK_IMPLEMENTATION" == "OpenBLAS" ]]; then
        pkgs="${pkgs} libopenblas-dev liblapacke-dev"
    elif [[ "$BLASLAPACK_IMPLEMENTATION" == "MKL" ]]; then
        # mkl_link_tool is a 32bits application !
        sudo dpkg --add-architecture i386
        sudo apt-get update
        pkgs="${pkgs} libc6:i386 libncurses5:i386 libstdc++6:i386"
    fi

    if [[ "$USE_VTK" == "1" ]]; then
        pkgs="${pkgs} libvtk5-dev"
    fi

    if [[ "$BUILD_DOCUMENTATION" == "ON" ]]; then
        pkgs="${pkgs} doxygen graphviz"
    fi

    if [[ "$USE_COVERAGE" == "1" ]]; then
        pkgs="${pkgs} lcov"
    fi

    if [[ "${pkgs}" != "" ]]; then
        sudo apt-get install -y ${pkgs}
    fi
fi

# install anaconda Python for wrapping or deployment
if [[ "$USE_PYTHON" == "1" ]]; then
    setup_conda_wrap
fi
