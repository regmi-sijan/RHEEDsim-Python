#!/bin/sh
# Portable launcher; prefers an installed Python with NumPy, Matplotlib and Tk.
cd "$(dirname "$0")" || exit 1
for candidate in ./.venv/bin/python python3 python3.13 python3.12 python3.11; do
    if "$candidate" -c 'import numpy, matplotlib, tkinter' >/dev/null 2>&1; then
        exec "$candidate" RHEEDsim.py "$@"
    fi
done
printf '%s\n' 'Install requirements.txt and a Python distribution with Tk, then run python3 RHEEDsim.py.'
exit 1
