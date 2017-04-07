if [[ "$USE_COVERAGE" == "1" ]]; then
    lcov --directory . --capture --output-file coverage.info # capture coverage info
    lcov --remove coverage.info '/usr/*' --output-file coverage.info # filter out system
    lcov --list coverage.info
    bash <(curl -s https://codecov.io/bash)
fi

