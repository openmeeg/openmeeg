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

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then

    # Install some custom requirements on OS X
    brew tap homebrew/science # a lot of cool formulae for scientific tools
    # brew update && brew upgrade

    # Install some custom requirements on OS X
    if [[ "$USE_PROJECT" == "0" || "$USE_SYSTEM" == "1" ]]; then
        brew install hdf5
        brew install libmatio
    fi

    if [[ "$BLASLAPACK_IMPLEMENTATION" == "OpenBLAS" || "$BLASLAPACK_IMPLEMENTATION" == "Auto" ]]; then
        brew install openblas
        brew link openblas --force  # required as link is not automatic
    fi

    if [[ "$USE_VTK" == "1" && "$ENABLE_PACKAGING" != "1" ]]; then
        brew install vtk
    fi

    if [[ "$USE_OMP" == 1 ]]; then
        brew install llvm
        export OMP_NUM_THREADS=4
        export LDFLAGS="-L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib"
        export DYLD_LIBRARY_PATH="/usr/local/opt/llvm/lib:$DYLD_LIBRARY_PATH"
        export CC="/usr/local/opt/llvm/bin/clang"
        export CXX="/usr/local/opt/llvm/bin/clang++"
        export CPPFLAGS="-I/usr/local/opt/llvm/include -fopenmp $CPPFLAGS"
        export CFLAGS="-I/usr/local/opt/llvm/include -fopenmp $CFLAGS"
    fi

    if [[ "$USE_CGAL" == 1 ]]; then
        brew install cgal
    fi

    if [[ "$BUILD_DOCUMENTATION" == "ON" ]]; then
        brew install Doxygen Graphviz # For building documentation
    fi

else
    # Install some custom requirements on Linux
    case "$COMPILER" in
    "gcc")
        export CXX=$(which g++)
        export CC=$(which gcc)
        ;;
    "clang")
        export CXX=$(which clang++)
        export CC=$(which clang)
        ;;
    "icc")
        wget -q -O install-icc.sh 'https://raw.githubusercontent.com/massich/icc-travis/install_only/install-icc.sh'
        chmod +x install-icc.sh
        ./install-icc.sh
        source ~/.bashrc
        echo "I'll install something here"
        echo "I'll install something and here"
        echo "I'll install something also"
        echo "I'll install something and more"
        export CXX=$(which icc)
        export CC=$(which icc)
        ;;
    *)
        echo "Unknown compiler"
        ;;
    esac

    echo "The compiler CXX is: ${CXX}"
    echo "The compiler CC is: ${CC}"

    if [[ "$USE_PROJECT" == "0" || "$USE_SYSTEM" == "1" ]]; then
        sudo apt-get install -y libhdf5-serial-dev libmatio-dev
    fi

    if [[ "$USE_CGAL" == 1 ]]; then
        sudo apt-get install -y libcgal-dev
    fi

    if [[ "$USE_GIFTI" == 1 ]]; then
        sudo apt-get install -y libnifti-dev libgiftiio-dev
    fi

    if [[ "$BLASLAPACK_IMPLEMENTATION" == "Atlas" ]]; then
        sudo apt-get install -y libatlas-dev libatlas-base-dev
    elif [[ "$BLASLAPACK_IMPLEMENTATION" == "LAPACK" ]]; then
        if [[ "$USE_PROJECT" == "0" || "$USE_SYSTEM" == "1" ]]; then
            sudo apt-get install -y liblapack-dev libblas-dev
        fi
    elif [[ "$BLASLAPACK_IMPLEMENTATION" == "OpenBLAS" ]]; then
        sudo apt-get install -y libopenblas-dev liblapacke-dev
    elif [[ "$BLASLAPACK_IMPLEMENTATION" == "MKL" ]]; then
        # mkl_link_tool is a 32bits application !
        echo "HERE"
        sudo dpkg --add-architecture i386
        sudo apt-get update
        sudo apt-get install -y libc6:i386 libncurses5:i386 libstdc++6:i386
    fi

    if [[ "$USE_VTK" == "1" ]]; then
        sudo apt-get install libvtk5-dev
    fi

    if [[ "$BUILD_DOCUMENTATION" == "ON" ]]; then
        sudo apt-get install -y doxygen graphviz
    fi

    if [[ "$USE_COVERAGE" == "1" ]]; then
        sudo apt-get install -y lcov
    fi
fi

# install anaconda Python for wrapping or deployment
if [[ "$USE_PYTHON" == "1" ]]; then
    setup_conda_wrap
fi
