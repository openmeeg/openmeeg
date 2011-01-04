#!/bin/bash

#PIPOL esn i386-linux-ubuntu-intrepid.dd.gz none 02:00 --user --silent
#PIPOL esn amd64-linux-ubuntu-intrepid.dd.gz none 02:00 --user --silent

#PIPOL esn i386-linux-ubuntu-jaunty.dd.gz none 02:00 --user --silent
#PIPOL esn amd64-linux-ubuntu-jaunty.dd.gz none 02:00 --user --silent

#PIPOL esn i386-linux-ubuntu-karmic.dd.gz none 02:00 --user --silent
#PIPOL esn amd64-linux-ubuntu-karmic.dd.gz none 02:00 --user --silent

#PIPOL esn i386-linux-fedora-core11.dd.gz none 02:00 --user --silent
#PIPOL esn amd64-linux-fedora-core11.dd.gz none 02:00 --user --silent

#PIPOL esn i386_kvm-linux-debian-lenny none 02:00 --user --silent
#PIPOL esn i386_kvm-linux-debian-testing none 02:00 --user --silent

#PIPOL esn amd64_kvm-linux-debian-lenny none 02:00 --user --silent
#PIPOL esn amd64_kvm-linux-debian-testing none 02:00 --user --silent

SYSTEM=`uname`

if [ -e ./pipol ] ; then
	rm -rf ./pipol/$PIPOL_HOST
	mkdir ./pipol/$PIPOL_HOST
else
	mkdir ./pipol
	rm -rf ./pipol/$PIPOL_HOST
	mkdir ./pipol/$PIPOL_HOST
fi
cd ./pipol/$PIPOL_HOST

svn checkout svn://scm.gforge.inria.fr/svn/openmeeg/trunk openmeeg-trunk --quiet

sh ./openmeeg-trunk/pipol/install_packages.sh
perl ./openmeeg-trunk/pipol/cmake.pl

cd openmeeg-trunk

rm CMakeCache.txt

#ucontext and pthread
cmake -DBUILD_TESTING=True -DENABLE_PACKAGING=True -DPYTHON_WRAP=True .
ctest -D NightlyStart
ctest -D NightlyConfigure
ctest -D NightlyBuild
ctest -D NightlyTest
ctest -D NightlySubmit
make clean
