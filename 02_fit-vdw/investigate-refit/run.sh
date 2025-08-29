#!/usr/bin/env bash
# source ~/.bashrc

# conda activate evaluator-050

# mkdir output images

# # collate data initially
# python collate-data.py -i ../refit -o output/training-per-iteration.csv

# # compare iterations from training for an initial look at improvement
# python compare-iterations-from-training.py -i output/training-per-iteration.csv -o images/iter_00-15.png
# python compare-iterations-from-training-N.py -i output/training-per-iteration.csv -o images/iter_00-15-N.png


# # compare parameters over iterations
# python compare-parameters-over-iterations.py -i ../refit -o images/parameters-over-time.png

# plot change in parameters
python plot-parameter-changes-percent.py                    \
    -if     openff-2.2.1.offxml                             \
    -of     ../refit/result/optimize/force-field.offxml     \
    -o      output                                          \
    -im     images
