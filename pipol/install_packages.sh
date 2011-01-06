#!/bin/bash

arch=`uname`

if [ -e /usr/bin/apt-get ] ; then
    sudo apt-get -y install subversion
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
    if [ x$PIPOL_IMAGE = xi386-linux-ubuntu-karmic.dd.gz ] ; then
    	sudo apt-get -y install valgrind
    fi
else
	if [ -e /usr/bin/yum ] ; then
		sudo yum -y update
	    sudo yum -y install subversion
	    sudo yum -y install gcc
	    sudo yum -y install make
	    sudo yum -y install cmake
	    sudo yum -y install wget
	    sudo yum -y install perl
        sudo yum -y install atlas-devel
        sudo yum -y install blas-devel
        sudo yum -y install numpy
        sudo yum -y install python-devel
	else
		if [ x$arch = xDarwin ] ; then
		    ###
		fi
	fi
fi

which_svn=`which svn`		#svn necessary
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
echo $which_svn
echo $which_wget
echo $which_perl
