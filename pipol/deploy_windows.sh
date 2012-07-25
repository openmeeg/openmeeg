#!/bin/bash

BRANCH=master
if [ x$1 == "x-release" ]; then
    BRANCH=$2
fi

scp openmeeg_windows.sh pipol:openmeeg_windows.sh

# win xp 32bits
ssh pipol pipol-sub i386_kvm-windows-xp-pro-sp3 03:00 "bash ~/openmeeg_windows.sh $BRANCH"

# win 7 64bits
ssh pipol pipol-sub amd64_kvm-windows-7 03:00 "bash ~/openmeeg_windows.sh $BRANCH"
