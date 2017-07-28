function setup_conda {
    if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
        curl https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh -o miniconda.sh -s
    else
        wget -q http://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh
    fi

    chmod +x miniconda.sh
    ./miniconda.sh -b -p ${HOME}/miniconda
    export PATH=${HOME}/miniconda/bin:$PATH
    conda update --yes --quiet conda
    export PYTHON=2.7
    conda create -n testenv --yes pip python=$PYTHON
    source activate testenv
}

if [[ "$USE_COVERAGE" == "1" ]]; then
    lcov --directory . --capture --output-file coverage.info > /dev/null 2>&1 # capture coverage info
    lcov --remove coverage.info '/usr/*' --output-file coverage.info > /dev/null 2>&1 # filter out system
    lcov --list coverage.info > /dev/null 2>&1
    bash <(curl -s https://codecov.io/bash) > /dev/null 2>&1
fi

# only upload to forge if we are on the master branch
if [[ $ENABLE_PACKAGING == "1" && $TRAVIS_PULL_REQUEST == "false" && $TRAVIS_BRANCH == "master" ]]; then
    if ! [ -x "$(command -v conda)" ]; then
        setup_conda
    fi
    conda install -y --quiet paramiko
    conda install -y --quiet pyopenssl
    python ${src_dir}/build_tools/upload_package_gforge.py *gz
fi
