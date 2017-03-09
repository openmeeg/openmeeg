function setup_conda {
    if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
        wget -q https://repo.continuum.io/miniconda/Miniconda-latest-MacOSX-x86_64.sh -O miniconda.sh
    else
        wget -q http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
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
    lcov --directory . --capture --output-file coverage.info # capture coverage info
    lcov --remove coverage.info '/usr/*' --output-file coverage.info # filter out system
    lcov --list coverage.info
    bash <(curl -s https://codecov.io/bash)
fi

if [[ $deploy_password != '' ]]; then
    echo $deploy_password > secret.txt
    openssl aes-256-cbc -pass "file:./secret.txt" -in $TRAVIS_BUILD_DIR/build_tools/travis/openmeeg_deploy_key.enc -out ./openmeeg_deploy_key -d -a
    chmod 600 openmeeg_deploy_key
    mv openmeeg_deploy_key ~/.ssh/openmeeg_deploy_key

    setup_conda
    conda install --yes --quiet paramiko
    python $TRAVIS_BUILD_DIR/build_tools/travis/upload_package_gforge.py *tar.gz
fi
