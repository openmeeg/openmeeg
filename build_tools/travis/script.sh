
"$CXX --version"

if [[ "$USE_PROJECT" == "0" ]]; then make CTEST_OUTPUT_ON_FAILURE=1 test; else make CTEST_OUTPUT_ON_FAILURE=1 test-OpenMEEG; fi

if [[ "$ENABLE_PACKAGING" == "1" ]]; then 
    cpack -G TGZ;

    # now test the binary package: uninstall all brew, extract previously built package, run om_assemble
    if [[ $APPLE_STANDALONE == "1" ]]; then
        /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall)"
        tar xvvzf OpenMEEG-2.*.gz
        cd OpenMEEG-2.*
        export DYLD_LIBRARY_PATH="lib:$DYLD_LIBRARY_PATH"
        ./bin/om_assemble
    fi
fi
