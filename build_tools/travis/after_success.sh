pwd

ls *

if [[ "$USE_COVERAGE" == "1" ]]; then
    bash <(curl -s https://codecov.io/bash)
fi
