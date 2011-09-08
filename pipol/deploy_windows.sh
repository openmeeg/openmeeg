#!/bin/bash

file=`mktemp -t tmp`
cp openmeeg_windows.sh $file
cat >> $file <<EOF
build_and_test "Experimental"
EOF

if [ x$1 == "x-release" ]; then
    # usage : ./deploy_windows.sh -release 2.1
    file2=`mktemp -t tmp`
    file3=`mktemp -t tmp`
    sed "s/openmeeg\/trunk/openmeeg\/branches\/release-$2/g" $file > $file2
    sed "s/openmeeg-trunk/openmeeg-release-$2/g" $file2 > $file3
    rm $file2
    mv $file3 $file
fi

chmod +x $file
scp $file pipol:experimental_windows.sh
rm $file

# win xp 32bits
ssh pipol pipol-sub esn i386_kvm-windows-xp-pro-sp3 none 01:00 "bash ~/experimental_windows.sh"

# win 7 64bits
ssh pipol pipol-sub esn amd64_kvm-windows-7 none 01:00 "bash ~/experimental_windows.sh"
