# The -e flag causes the script to exit as soon as one command returns a non-zero exit code.
set -e

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    find . -name "om_assemble" | xargs otool -L
    find . -name "libOpenMEEG*dylib" | xargs otool -L
else
    find . -name "om_assemble" | xargs ldd
    find . -name "libOpenMEEG*so" | xargs ldd
fi

if [[ "$USE_PROJECT" == "0" ]]; then
    make CTEST_OUTPUT_ON_FAILURE=1 test;
else
    make CTEST_OUTPUT_ON_FAILURE=1 test-OpenMEEG;
fi

if [[ "$ENABLE_PACKAGING" == "1" ]]; then
    cpack -G TGZ;

    if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
        # remove completly brew to test standalone mac package
        curl -O https://raw.githubusercontent.com/Homebrew/install/master/uninstall
        chmod +x uninstall
        ./uninstall --force
    else
        sudo apt-get remove -y libhdf5-serial-dev libmatio-dev libopenblas-dev liblapacke-dev libvtk5-dev libvtk5.8 libnifti-dev libgiftiio-dev libcgal-dev doxygen graphviz
    fi

    # Remove conda used from wrapping
    echo "Removing miniconda_wrap"
    rm -rf ${HOME}/miniconda_wrap

    if [[ "$BLASLAPACK_IMPLEMENTATION" == "MKL" ]]; then
        echo "Removing MKL"
        sudo rm -rf /opt/intel
    fi

    echo "Untaring the package"
    tar xzf OpenMEEG-2.*.gz #> /dev/null 2> /dev/null
    cd OpenMEEG-2.*
    if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
        echo "Running otool on om_assemble"
        otool -L ./bin/om_assemble
    else
        echo "Running ldd on om_assemble"
        ldd ./bin/om_assemble
    fi
    echo "Setting up (DY)LD_LIBRARY_PATH"
    export DYLD_LIBRARY_PATH="lib:$DYLD_LIBRARY_PATH"
    export LD_LIBRARY_PATH="lib:$LD_LIBRARY_PATH"
    echo "running ./bin/om_assemble"
    ./bin/om_assemble
    cd ..
fi
