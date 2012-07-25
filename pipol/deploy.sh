#!/bin/bash

# usage : ./deploy.sh               # for a simple experimental build
#         ./deploy.sh -full         # for a full experimental build
#         ./deploy.sh -release 2.1  # for a release

BRANCH=master
if [ x$1 == "x-release" ]; then
    BRANCH=$2
fi

scp openmeeg_unix.sh pipol:openmeeg_unix.sh
scp install_packages.sh pipol:install_packages.sh

# linux 32
ssh pipol pipol-sub i386_kvm-linux-debian-lenny 02:00 "~/openmeeg_unix.sh $BRANCH"
ssh pipol pipol-sub i386-linux-fedora-core11.dd.gz 02:00 "~/openmeeg_unix.sh $BRANCH"

# linux 64
ssh pipol pipol-sub amd64_kvm-linux-debian-lenny 02:00 "~/openmeeg_unix.sh $BRANCH"
ssh pipol pipol-sub amd64-linux-fedora-core11.dd.gz 02:00 "~/openmeeg_unix.sh $BRANCH"

# mac 32
ssh pipol pipol-sub i386_mac-mac-osx-server-leopard 02:00 "~/openmeeg_unix.sh $BRANCH"

# mac 64
ssh pipol pipol-sub x86_64_mac-mac-osx-server-snow-leopard 02:00 "~/openmeeg_unix.sh $BRANCH"

if [ x$1 == "x-full" ]; then
    ssh pipol pipol-sub linux 02:00 "~/openmeeg_unix.sh"
    ssh pipol pipol-sub mac 02:00 "~/openmeeg_unix.sh"
    # ssh pipol pipol-sub i386-unix-solaris-10 02:00 "~/openmeeg_unix.sh"
fi
