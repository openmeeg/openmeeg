if [[ $TRAVIS_OS_NAME == 'osx' ]]; then

    # Install some custom requirements on OS X
    brew tap homebrew/science # a lot of cool formulae for scientific tools
    brew tap homebrew/python # numpy, scipy, matplotlib, ...
    brew update && brew upgrade


    # Install some custom requirements on OS X
    if [[ "$USE_PROJECT" == "0" ]]; then
        brew install hdf5
    fi

    # install a brewed python
    # brew install python
    # brew install numpy

    brew install cmake
    brew install swig
    brew install Doxygen

    if [[ "$USE_OPENBLAS" == "1" ]]; then
        brew install openblas
    fi

else
    # Install some custom requirements on Linux
    # g++4.8.1
    if [ "$CXX" == "g++" ]; then sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test; fi

    # clang 3.4
    if [ "$CXX" == "clang++" ]; then sudo add-apt-repository -y ppa:h-rayflood/llvm; fi

    if [[ "$USE_SYSTEM" == "1" ]]; then
        # To get recent hdf5 version:
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
    fi

    sudo apt-get update -qq
fi
