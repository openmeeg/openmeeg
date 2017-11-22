if [[ "$USE_COVERAGE" == "1" ]]; then
    lcov --directory . --capture --output-file coverage.info > /dev/null 2>&1 # capture coverage info
    lcov --remove coverage.info '/usr/*' --output-file coverage.info > /dev/null 2>&1 # filter out system
    lcov --list coverage.info > /dev/null 2>&1
    bash <(curl -s https://codecov.io/bash) > /dev/null 2>&1
fi


# only upload to forge if we are on the master branch
if [[ $ENABLE_PACKAGING == "1" && $TRAVIS_PULL_REQUEST == "false" && $TRAVIS_BRANCH == "master" ]]; then
    if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
        curl https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh -o miniconda.sh -s
    else
        wget -q http://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh
    fi

    chmod +x miniconda.sh
    ./miniconda.sh -b -p ${HOME}/miniconda_deploy

    export PATH=${HOME}/miniconda_deploy/bin:$PATH
    conda update --yes --quiet conda
    conda create -n packageenv --yes pip python=2.7
    source activate packageenv
    conda install -y --quiet paramiko
    conda install -y --quiet pyopenssl

    python ${src_dir}/build_tools/upload_package_gforge.py *gz
fi

# tear down ICC compiler

if [[ "$COMPILER" == "icc" ]]; then
    echo "[SK-DB-MSG] TEARING DOWN ICC"
    [[ ! -z "${INTEL_INSTALL_PATH}" ]] && uninstall_intel_software
fi
