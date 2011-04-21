#!/bin/bash

scp experimental.sh pipol:.

# linux 32
ssh pipol pipol-sub esn i386_kvm-linux-debian-lenny none 02:00 "~/experimental.sh"
ssh pipol pipol-sub esn i386-linux-fedora-core11.dd.gz none 02:00 "~/experimental.sh"
# ssh pipol pipol-sub esn i386_kvm-linux-debian-testing none 02:00 "~/experimental.sh"

# linux 64
ssh pipol pipol-sub esn amd64_kvm-linux-debian-lenny none 02:00 "~/experimental.sh"
ssh pipol pipol-sub esn amd64-linux-fedora-core11.dd.gz none 02:00 "~/experimental.sh"
# ssh pipol pipol-sub esn amd64_kvm-linux-debian-testing none 02:00 "~/experimental.sh"

# # mac 32
ssh pipol pipol-sub esn i386_mac-mac-osx-server-leopard none 02:00 "~/experimental.sh"
# mac 64
# ssh pipol pipol-sub esn x86_64_mac-mac-osx-server-snow-leopard none 02:00 "~/experimental.sh"
