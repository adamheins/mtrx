#!/bin/bash

# Run clar mixer to generate test metadata.
python ../../clar/generate.py -f .

# Compile the project and tests.
gcc clar.c main.c mtrx_suite.c ../mtrx.c -lm -o testit -std=c99 -D_GNU_SOURCE

# Run the tests.
./testit

# Remove the test binary.
rm testit