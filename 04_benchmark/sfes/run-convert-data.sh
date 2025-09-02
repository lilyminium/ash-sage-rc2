#!/bin/bash

python convert-mnsol-data.py > logs/convert-data.log
python convert-freesolv-data.py >> logs/convert-data.log