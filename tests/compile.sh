#!/bin/bash
gcc clar.c main.c mtrx_suite.c ../mtrx.c -lm -o testit -std=c99 -D_GNU_SOURCE
