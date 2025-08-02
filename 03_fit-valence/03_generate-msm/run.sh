
#!/bin/bash


FFDIR="../01_generate-forcefield/output"
DATASET="../../01_download-data/qm/hessian-data/hessian_results.json"


mkdir -p output logs msm-ff

# conda activate sage-2.3.0-qubekit
micromamba env export > full-qubekit-env.yaml

python calculate-msm.py                         \
    -i      $DATASET                            \
    -o      msm-data > logs/calculate-msm.log 2>&1


for ver in v0 v1 v2 v3 ; do

    echo "==== Generating FF for ${ver} ===="
    python generate-msm-forcefield.py       \
        -im     msm-data                    \
        -om     "msm-ff/ff-${ver}.json"     \
        -i      "${FFDIR}/initial-force-field-${ver}.offxml"    \
        -a      mean                                        \
        -o      output/msm-force-field-${ver}-mean.offxml > msm-${ver}-mean.log 2>&1

done

