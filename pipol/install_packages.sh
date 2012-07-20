#!/bin/bash

arch=`uname`

if [ -e /usr/bin/apt-get ] ; then
    sudo apt-get -y install git-core
    sudo apt-get -y install gcc
    sudo apt-get -y install g++
    sudo apt-get -y install make
    sudo apt-get -y install cmake
    sudo apt-get -y install wget
    sudo apt-get -y install perl
    sudo apt-get -y install libatlas-base-dev
    sudo apt-get -y install python-numpy
    sudo apt-get -y install python-dev
    sudo apt-get -y install swig

    # For git
    sudo apt-get -y install libssl-dev
    sudo apt-get -y install libcurl4-openssl-dev
    sudo apt-get -y install libexpat1-dev
    sudo apt-get -y install perl-modules
    sudo apt-get -y install gettext

    # For memory checking

    if [ x$PIPOL_IMAGE = xi386-linux-ubuntu-karmic.dd.gz ] ; then
    	sudo apt-get -y install valgrind
    fi
else
	if [ -e /usr/bin/yum ] ; then
		sudo yum -y update
	    sudo yum -y install git-core
	    sudo yum -y install gcc
	    sudo yum -y install make
	    sudo yum -y install cmake
	    sudo yum -y install wget
	    sudo yum -y install perl
        sudo yum -y install atlas-devel
        sudo yum -y install blas-devel
        sudo yum -y install numpy
        sudo yum -y install swig
        sudo yum -y install python-devel

        # For git
        sudo yum -y install openssl-devel
        sudo yum -y install libcurl-devel
        sudo yum -y install expat-devel
        sudo yum -y install perl-ExtUtils-MakeMaker
    else
       if [ x$arch = xDarwin ] ; then
            /usr/bin/ruby -e "$(/usr/bin/curl -fsSL https://raw.github.com/mxcl/homebrew/master/Library/Contributions/install_homebrew.rb)"
            brew install git
       fi
    fi
fi

# Old versions of git (before 1.6.5) do not know the --recursive option...

my_git_version=`git --version`
case $my_git_version in
    *1.[7-9].[0-9]*)
        ;;
    *1.6.[5-9]*)
        ;;
    *)
        wget http://git-core.googlecode.com/files/git-1.7.11.2.tar.gz
        tar zxvf git-1.7.11.2.tar.gz
        cd ./git-1.7.11.2
        make prefix=/usr all
        if [ -e /usr/bin/apt-get ] ; then
            sudo apt-get -y remove git-core
        else
            if [ -e /usr/bin/yum ] ; then
                sudo yum -y remove git-core
            fi
        fi
        sudo make prefix=/usr install;;
esac

which_git=`which git`		#git necessary
which_gcc=`which gcc`		#gcc gcc necessary
which_gpp=`which g++`		#gcc g++ necessary
which_make=`which make`		#make necessary
which_cmake=`which cmake`	#cmake necessary
which_wget=`which wget`		#wget for cmake
which_perl=`which perl`		#perl
echo $which_cmake
echo $which_unzip
echo $which_make
echo $which_gcc
echo $which_gpp
echo $which_git
echo $which_wget
echo $which_perl
