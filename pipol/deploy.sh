#!/bin/bash

scp experimental.sh pipol:.

# linux 32
ssh pipol pipol-sub esn i386_kvm-linux-debian-lenny none 02:00 "~/experimental.sh"
ssh pipol pipol-sub esn i386-linux-fedora-core11.dd.gz none 02:00 "~/experimental.sh"

# linux 64
ssh pipol pipol-sub esn amd64_kvm-linux-debian-lenny none 02:00 "~/experimental.sh"
ssh pipol pipol-sub esn amd64-linux-fedora-core11.dd.gz none 02:00 "~/experimental.sh"

# mac 32
ssh pipol pipol-sub esn i386_mac-mac-osx-server-leopard none 02:00 "~/experimental.sh"

# mac 64
#ssh pipol pipol-sub esn x86_64_mac-mac-osx-server-snow-leopard none 02:00 "~/experimental.sh"

if [ x$1 != "x-simple" ]; then

    ssh pipol pipol-sub esn i386-linux-ubuntu-intrepid.dd.gz none 02:00 "~/experimental.sh"
    ssh pipol pipol-sub esn amd64-linux-ubuntu-intrepid.dd.gz none 02:00 "~/experimental.sh"

    ssh pipol pipol-sub esn i386-linux-ubuntu-jaunty.dd.gz none 02:00 "~/experimental.sh"
    ssh pipol pipol-sub esn amd64-linux-ubuntu-jaunty.dd.gz none 02:00 "~/experimental.sh"

    ssh pipol pipol-sub esn i386-linux-ubuntu-karmic.dd.gz none 02:00 "~/experimental.sh"
    ssh pipol pipol-sub esn amd64-linux-ubuntu-karmic.dd.gz none 02:00 "~/experimental.sh"

    ssh pipol pipol-sub esn i386-linux-ubuntu-lucid.dd.gz none 02:00 "~/experimental.sh"
    ssh pipol pipol-sub esn amd64-linux-ubuntu-lucid.dd.gz none 02:00 "~/experimental.sh"

    ssh pipol pipol-sub esn i386_kvm-linux-debian-testing none 02:00 "~/experimental.sh"
    ssh pipol pipol-sub esn amd64_kvm-linux-debian-testing none 02:00 "~/experimental.sh"

    ssh pipol pipol-sub esn i386-linux-fedora-core12.dd.gz none 02:00 "~/experimental.sh"
    ssh pipol pipol-sub esn i386-linux-fedora-core13.dd.gz none 02:00 "~/experimental.sh"
    ssh pipol pipol-sub esn i386_2010-linux-fedora-core14.dd.gz none 02:00 "~/experimental.sh"

    ssh pipol pipol-sub esn amd64-linux-fedora-core12.dd.gz none 02:00 "~/experimental.sh"
    ssh pipol pipol-sub esn amd64-linux-fedora-core13.dd.gz none 02:00 "~/experimental.sh"
    ssh pipol pipol-sub esn amd64_2010-linux-fedora-core14.dd.gz none 02:00 "~/experimental.sh"

    ssh pipol pipol-sub esn amd64-linux-opensuse-11.dd.gz none 02:00 "~/experimental.sh"

    # ssh pipol pipol-sub esn i386-unix-solaris-10 none 02:00 "~/experimental.sh"

    ssh pipol pipol-sub esn i386-linux-centos-5 none 02:00 "~/experimental.sh"
    ssh pipol pipol-sub esn amd64-linux-centos-5 none 02:00 "~/experimental.sh"
fi
