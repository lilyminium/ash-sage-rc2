#!/bin/bash

conda activate evaluator-050

mkdir logs

# python remap-from-rc1-results.py \
#     -old ../../../old-ash-sage/01_download-data/physprop/final/output/validation-set.json \
#     -new ../../01_download-data/physprop/final/output/validation-set.json \
#     -o mappings/old-validation-to-new-validation.json > logs/remap-validation-to-validation.log

# python remap-from-rc1-results.py \
#     -old ../../../old-ash-sage/01_download-data/physprop/final/output/validation-set.json \
#     -new ../../01_download-data/physprop/final/output/training-set.json \
#     -o mappings/old-validation-to-new-training.json > logs/remap-validation-to-training.log

# python remap-from-rc1-results.py \
#     -old ../../../old-ash-sage/01_download-data/physprop/final/output/training-set.json \
#     -new ../../01_download-data/physprop/final/output/validation-set.json \
#     -o mappings/old-training-to-new-validation.json > logs/remap-training-to-validation.log

# python remap-from-rc1-results.py \
#     -old ../../../old-ash-sage/01_download-data/physprop/final/output/training-set.json \
#     -new ../../01_download-data/physprop/final/output/training-set.json \
#     -o mappings/old-training-to-new-training.json > logs/remap-training-to-training.log

RC1_DIRECTORY="../../../ash-sage/04_benchmark/phys-prop/thermoml-validation"

python copy-files-from-rc1.py -i "${RC1_DIRECTORY}/validation" -o validation mappings/old-validation-to-new-validation.json
python copy-files-from-rc1.py -i "${RC1_DIRECTORY}/training" -o validation mappings/old-training-to-new-validation.json
python copy-files-from-rc1.py -i "${RC1_DIRECTORY}/validation" -o training mappings/old-validation-to-new-training.json
python copy-files-from-rc1.py -i "${RC1_DIRECTORY}/training" -o training mappings/old-training-to-new-training.json
