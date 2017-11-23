# The -e flag causes the script to exit as soon as one command returns a non-zero exit code.
set -e

if [[ "$USE_PROJECT" == "0" ]]; then
    make CTEST_OUTPUT_ON_FAILURE=1 test;
else
    make CTEST_OUTPUT_ON_FAILURE=1 test-OpenMEEG;
fi

if [[ "$ENABLE_PACKAGING" == "1" ]]; then
    cpack -G TGZ;

    if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
        # remove completly brew to test standalone mac package
        /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall)"
    else
        sudo apt-get remove -y libhdf5-serial-dev libmatio-dev libopenblas-dev liblapacke-dev libvtk5-dev libvtk5.8
    fi

    # Remove conda used from wrapping
    rm -rf ${HOME}/miniconda_wrap

    if [[ "$BLASLAPACK_IMPLEMENTATION" == "MKL" ]]; then
        sudo rm -rf /opt/intel
    fi

    tar xvvzf OpenMEEG-2.*.gz > /dev/null 2> /dev/null
    cd OpenMEEG-2.*
    export DYLD_LIBRARY_PATH="lib:$DYLD_LIBRARY_PATH"
    export LD_LIBRARY_PATH="lib:$LD_LIBRARY_PATH"
    echo "running ./bin/om_assemble"
    ./bin/om_assemble
    cd ..
fi
