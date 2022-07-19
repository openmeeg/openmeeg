#!/bin/bash

set -xe
pytest $(python -c 'from pathlib import Path; import openmeeg; print(Path(openmeeg.__file__).parent)')