#!/usr/bin/env bash

python generate-comparison-patterns-from-ff.py \
    -i ../forcefields/fb-fit-v3-single-mean-k100.offxml \
    -o comparison-patterns/patterns-v3.json > logs/generate-comparison-patterns-v3.log
