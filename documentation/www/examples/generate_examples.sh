#!/usr/bin/env bash

pygmentize -O full,style=default -o python.html compute_leadfields.py
pygmentize -O full,style=default -o bash_script.html compute_leadfields.sh
pygmentize -O full,style=default -o windows_bat.html compute_leadfields.bat
