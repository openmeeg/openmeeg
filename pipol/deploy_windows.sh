#!/bin/bash

scp openmeeg_windows.sh pipol:.

# win xp 32bits
ssh pipol pipol-sub esn i386_kvm-windows-xp-pro-sp3 none 00:40 "bash ~/openmeeg_windows.sh Experimental"

# win 7 64bits
ssh pipol pipol-sub esn amd64_kvm-windows-7 none 00:40 "bash ~/openmeeg_windows.sh Experimental"
