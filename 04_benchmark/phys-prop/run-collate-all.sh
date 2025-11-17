#!/bin/bash

for tier in training validation ; do
    python collate-all-results.py   \
    -i results/${tier}              \
    -o output/${tier} \
    > logs/collate/collate-phys-prop-${tier}.log
done

