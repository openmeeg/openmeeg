if [[ $TRAVIS_OS_NAME == 'osx' ]]; then

    # Install some custom requirements on OS X
    brew tap homebrew/science # a lot of cool formulae for scientific tools
    brew tap homebrew/python # numpy, scipy, matplotlib, ...
    brew update && brew upgrade

    # Install some custom requirements on OS X
    if [[ "$USE_PROJECT" == "0" || "$USE_SYSTEM" == "1" ]]; then
        brew install hdf5
        brew install libmatio --with-hdf5
    fi

    # install a brewed python
    # To use Python of
    if [[ "$USE_PYTHON" == "1" ]]; then
        brew install python
        brew install numpy
        brew install swig
    fi

    # brew install Doxygen  # For building documentation

    if [[ "$USE_OPENBLAS" == "1" ]]; then
        # brew install liblapacke ?
        brew install openblas
        brew link openblas --force  # required as link is not automatic
    fi

    if [[ "$USE_VTK" == "1" ]]; then
        brew install vtk
    fi

    brew install cmake

else
    # Install some custom requirements on Linux
    # g++4.8.1
    if [ "$CXX" == "g++" ]; then sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test; fi

    # clang 3.4
    if [ "$CXX" == "clang++" ]; then sudo add-apt-repository -y ppa:h-rayflood/llvm; fi

    # To get recent hdf5 version: AND openblas
    sudo sed -i -e 's/trusty/vivid/g' /etc/apt/sources.list
        # # Handle MATIO
        # sudo apt-get update -qq
        # # to prevent IPv6 being used for APT
        # sudo bash -c "echo 'Acquire::ForceIPv4 \"true\";' > /etc/apt/apt.conf.d/99force-ipv4"
        # # The ultimate one-liner setup for NeuroDebian repository
        # bash <(wget -q -O- http://neuro.debian.net/_files/neurodebian-travis.sh)
        # # But we actually want -devel repository just to get matio backport
        # sudo sed -ie 's,neuro.debian.net/debian ,neuro.debian.net/debian-devel ,g' /etc/apt/sources.list.d/neurodebian.sources.list
        # # Just to get information about available versions
        # apt-cache policy libmatio-dev

    sudo apt-get update -qq

    # # g++4.8.1
    if [ "$CXX" = "g++" ]; then sudo apt-get install -qq g++-4.8; fi
    if [ "$CXX" = "g++" ]; then export CXX="g++-4.8"; fi

    # clang 3.4
    if [ "$CXX" == "clang++" ]; then sudo apt-get install --allow-unauthenticated -qq clang-3.4; fi
    if [ "$CXX" == "clang++" ]; then export CXX="clang++-3.4"; fi

    sudo apt-get update
    if [[ "$USE_SYSTEM" == "1" ]]; then
      sudo apt-get install -y libhdf5-serial-dev libopenblas-base
    fi
    if [[ "$USE_ATLAS" == "1" ]]; then
        sudo apt-get install -y libatlas-dev libatlas-base-dev libblas-dev liblapack-dev
    fi
    if [[ "$USE_OPENBLAS" == "1" ]]; then
        sudo apt-get install -y libopenblas-dev
    fi

    if [[ "$USE_PYTHON" == "1" ]]; then
        sudo apt-get install -y swig python-dev python-numpy
    fi

    if [[ "$USE_VTK" == "1" ]]; then
        sudo apt-get install libvtk5-dev libtiff4-dev
    fi

    # sudo apt-get install doxygen
    wget https://s3.amazonaws.com/biibinaries/thirdparty/cmake-3.0.2-Linux-64.tar.gz
    tar -xzf cmake-3.0.2-Linux-64.tar.gz
    sudo cp -fR cmake-3.0.2-Linux-64/* /usr
    rm -rf cmake-3.0.2-Linux-64
    rm cmake-3.0.2-Linux-64.tar.gz
fi
