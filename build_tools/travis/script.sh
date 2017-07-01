
"$CXX --version"

if [[ "$USE_PROJECT" == "0" ]]; then
    make CTEST_OUTPUT_ON_FAILURE=1 test;
else
    make CTEST_OUTPUT_ON_FAILURE=1 test-OpenMEEG;
fi

if [[ "$ENABLE_PACKAGING" == "1" ]]; then
    cpack -G TGZ;

    # now test the binary package: uninstall all brew, extract previously built package, run om_assemble
    if [[ $APPLE_STANDALONE == "1" ]]; then
        # do not completly remove brew we need wget, ... for conda
        # /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall)"
        brew uninstall --ignore-dependencies --force hdf5 libmatio vtk cgal openblas
        tar xvvzf OpenMEEG-2.*.gz > /dev/null 2> /dev/null
        cd OpenMEEG-2.*
        export DYLD_LIBRARY_PATH="lib:$DYLD_LIBRARY_PATH"
        ./bin/om_assemble
        cd ..
    fi
fi
