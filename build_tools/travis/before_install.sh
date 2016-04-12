if [[ $TRAVIS_OS_NAME == 'osx' ]]; then

    # Install some custom requirements on OS X
    brew tap homebrew/science # a lot of cool formulae for scientific tools
    brew tap homebrew/python # numpy, scipy, matplotlib, ...
    brew update && brew upgrade

    # install a brewed python
    brew install python
    brew install numpy
    brew install cmake
    brew install swig
    brew install Doxygen
    brew install https://raw.githubusercontent.com/kadwanev/bigboybrew/master/Library/Formula/sshpass.rb

else
    # Install some custom requirements on Linux
    # g++4.8.1
    if [ "$CXX" == "g++" ]; then sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test; fi

    # clang 3.4
    if [ "$CXX" == "clang++" ]; then sudo add-apt-repository -y ppa:h-rayflood/llvm; fi

    sudo apt-get update -qq

fi
