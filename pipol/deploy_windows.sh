#!/bin/bash

scp experimental_windows.sh pipol:.

# win xp 32bits
ssh pipol pipol-sub esn i386_kvm-windows-xp-pro-sp3 none 01:00 "bash ~/experimental_windows.sh"

# win 7 64bits
ssh pipol pipol-sub esn amd64_kvm-windows-7 none 01:00 "bash ~/experimental_windows.sh"
