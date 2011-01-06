#!/bin/bash

scp experimental.sh pipol:.

ssh pipol pipol-sub esn i386_kvm-linux-debian-lenny none 02:00 "~/experimental.sh"
ssh pipol pipol-sub esn i386-linux-fedora-core11.dd.gz none 02:00 "~/experimental.sh"
# ssh pipol pipol-sub esn i386_kvm-linux-debian-testing none 02:00 "~/experimental.sh"

ssh pipol pipol-sub esn amd64_kvm-linux-debian-lenny none 02:00 "~/experimental.sh"
ssh pipol pipol-sub esn amd64-linux-fedora-core11.dd.gz none 02:00 "~/experimental.sh"
# ssh pipol pipol-sub esn amd64_kvm-linux-debian-testing none 02:00 "~/experimental.sh"
