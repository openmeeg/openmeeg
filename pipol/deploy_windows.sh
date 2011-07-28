#!/bin/bash

file=`mktemp`
cp openmeeg_windows.sh $file
cat >> $file <<EOF
build_and_test "Experimental"
EOF

chmod +x $file
scp $file pipol:experimental_windows.sh
rm $file

# win xp 32bits
ssh pipol pipol-sub esn i386_kvm-windows-xp-pro-sp3 none 01:00 "bash ~/experimental_windows.sh"

# win 7 64bits
ssh pipol pipol-sub esn amd64_kvm-windows-7 none 01:00 "bash ~/experimental_windows.sh"
