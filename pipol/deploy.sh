#!/bin/bash

# usage : ./deploy.sh               # for a simple experimental build
#         ./deploy.sh -full         # for a full experimental build
#         ./deploy.sh -release 2.1  # for a release


file=`mktemp`
cp openmeeg_unix.sh $file

if [ x$1 == "x-release" ]; then
    file2=`mktemp`
    file3=`mktemp`
    sed "s/openmeeg\/trunk/openmeeg\/branches\/release-$2/g" $file > $file2
    sed "s/openmeeg-trunk/openmeeg-release-$2/g" $file2 > $file3
    rm $file2
    mv $file3 $file
fi

chmod +x $file
scp $file pipol:openmeeg_unix.sh
rm $file

# linux 32
ssh pipol pipol-sub esn i386_kvm-linux-debian-lenny none 02:00 "~/openmeeg_unix.sh"
ssh pipol pipol-sub esn i386-linux-fedora-core11.dd.gz none 02:00 "~/openmeeg_unix.sh"

# linux 64
ssh pipol pipol-sub esn amd64_kvm-linux-debian-lenny none 02:00 "~/openmeeg_unix.sh"
ssh pipol pipol-sub esn amd64-linux-fedora-core11.dd.gz none 02:00 "~/openmeeg_unix.sh"

# mac 32
ssh pipol pipol-sub esn i386_mac-mac-osx-server-leopard none 02:00 "~/openmeeg_unix.sh"

# mac 64
ssh pipol pipol-sub esn x86_64_mac-mac-osx-server-snow-leopard none 02:00 "~/openmeeg_unix.sh"


if [ x$1 == "x-full" ]; then

    ssh pipol pipol-sub esn i386-linux-ubuntu-intrepid.dd.gz none 02:00 "~/openmeeg_unix.sh"
    ssh pipol pipol-sub esn amd64-linux-ubuntu-intrepid.dd.gz none 02:00 "~/openmeeg_unix.sh"

    ssh pipol pipol-sub esn i386-linux-ubuntu-jaunty.dd.gz none 02:00 "~/openmeeg_unix.sh"
    ssh pipol pipol-sub esn amd64-linux-ubuntu-jaunty.dd.gz none 02:00 "~/openmeeg_unix.sh"

    ssh pipol pipol-sub esn i386-linux-ubuntu-karmic.dd.gz none 02:00 "~/openmeeg_unix.sh"
    ssh pipol pipol-sub esn amd64-linux-ubuntu-karmic.dd.gz none 02:00 "~/openmeeg_unix.sh"

    ssh pipol pipol-sub esn i386-linux-ubuntu-lucid.dd.gz none 02:00 "~/openmeeg_unix.sh"
    ssh pipol pipol-sub esn amd64-linux-ubuntu-lucid.dd.gz none 02:00 "~/openmeeg_unix.sh"

    ssh pipol pipol-sub esn i386_kvm-linux-debian-testing none 02:00 "~/openmeeg_unix.sh"
    ssh pipol pipol-sub esn amd64_kvm-linux-debian-testing none 02:00 "~/openmeeg_unix.sh"

    ssh pipol pipol-sub esn i386-linux-fedora-core12.dd.gz none 02:00 "~/openmeeg_unix.sh"
    ssh pipol pipol-sub esn i386-linux-fedora-core13.dd.gz none 02:00 "~/openmeeg_unix.sh"
    ssh pipol pipol-sub esn i386_2010-linux-fedora-core14.dd.gz none 02:00 "~/openmeeg_unix.sh"

    ssh pipol pipol-sub esn amd64-linux-fedora-core12.dd.gz none 02:00 "~/openmeeg_unix.sh"
    ssh pipol pipol-sub esn amd64-linux-fedora-core13.dd.gz none 02:00 "~/openmeeg_unix.sh"
    ssh pipol pipol-sub esn amd64_2010-linux-fedora-core14.dd.gz none 02:00 "~/openmeeg_unix.sh"

    ssh pipol pipol-sub esn amd64-linux-opensuse-11.dd.gz none 02:00 "~/openmeeg_unix.sh"

    # ssh pipol pipol-sub esn i386-unix-solaris-10 none 02:00 "~/openmeeg_unix.sh"

    ssh pipol pipol-sub esn i386-linux-centos-5 none 02:00 "~/openmeeg_unix.sh"
    ssh pipol pipol-sub esn amd64-linux-centos-5 none 02:00 "~/openmeeg_unix.sh"
fi
