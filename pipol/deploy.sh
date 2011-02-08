#!/bin/bash

scp experimental.sh pipol:.

ssh pipol pipol-sub esn i386-linux-ubuntu-intrepid.dd.gz none 02:00 "~/experimental.sh"
ssh pipol pipol-sub esn amd64-linux-ubuntu-intrepid.dd.gz none 02:00 "~/experimental.sh"

ssh pipol pipol-sub esn i386-linux-ubuntu-jaunty.dd.gz none 02:00 "~/experimental.sh"
ssh pipol pipol-sub esn amd64-linux-ubuntu-jaunty.dd.gz none 02:00 "~/experimental.sh"

ssh pipol pipol-sub esn i386-linux-ubuntu-karmic.dd.gz none 02:00 "~/experimental.sh"
ssh pipol pipol-sub esn amd64-linux-ubuntu-karmic.dd.gz none 02:00 "~/experimental.sh"

ssh pipol pipol-sub esn i386-linux-fedora-core11.dd.gz none 02:00 "~/experimental.sh"
ssh pipol pipol-sub esn amd64-linux-fedora-core11.dd.gz none 02:00 "~/experimental.sh"

ssh pipol pipol-sub esn i386_kvm-linux-debian-lenny none 02:00 "~/experimental.sh"
ssh pipol pipol-sub esn i386_kvm-linux-debian-testing none 02:00 "~/experimental.sh"

ssh pipol pipol-sub esn amd64_kvm-linux-debian-lenny none 02:00 "~/experimental.sh"
ssh pipol pipol-sub esn amd64_kvm-linux-debian-testing none 02:00 "~/experimental.sh"

ssh pipol pipol-sub esn x86_64_mac-mac-osx-server-snow-leopard none 02:00 "~/experimental.sh"
# ssh pipol pipol-sub esn x86_mac-mac-osx-server-snow-leopard none 02:00 "~/experimental.sh"

ssh pipol pipol-sub esn i386-linux-centos-5 none 02:00 "~/experimental.sh"
ssh pipol pipol-sub esn amd64-linux-centos-5 none 02:00 "~/experimental.sh"
