#!/bin/sh
# Reruns the restriction-step certificates of Theorem B(i). Needs python3 and mpmath; about a minute.
set -e
cd "$(dirname "$0")"
python3 small_dimensions.py --max-n 15 --output small_certificate.json > small_dimensions.out
python3 large_cutoff_certificate.py > large_cutoff_certificate.out
python3 certificate.py > certificate.out
python3 audit.py > audit.out
printf 'all restriction certificates passed\n'
