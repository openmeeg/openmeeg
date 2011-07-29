#!/bin/bash

# script to update the openmeeg nightly build config on pipol

file=`mktemp`
cp openmeeg_windows.sh $file
cat >> $file <<EOF
build_and_test "Nightly"
EOF

chmod +x $file
scp $file pipol:~/.pipol/nightly/openmeeg_windows_nightly.sh
rm $file

scp openmeeg_nightly.sh pipol:~/.pipol/nightly/
