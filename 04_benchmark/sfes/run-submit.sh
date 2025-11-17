#!/bin/bash

# conda activate pontibus-alchemiscale-022

mkdir scoped-keys

# FFNAME="fb-fit-v3-single-mean-k100"
# FFNAME="fb-fit-v3-single-mean-k20"
# FFNAME="fit-iter-1"
FFNAME="vdw-refit"


FFNAME_UNDERSCORED=${FFNAME//-/_}

echo $FFNAME_UNDERSCORED

python submit.py                                                    \
    --org_scope "openff" --repeats 3                                \
    --network_filename      networks/mnsol-${FFNAME}-network.json   \
    --scope_name_campaign   "${FFNAME_UNDERSCORED}"                 \
    --scopekey_output       scoped-keys/mnsol_${FFNAME}_key.dat     \
    --scope_name_project    "mnsol" > logs/submit/submit-mnsol-${FFNAME}.log

python submit.py                                                    \
    --org_scope "openff" --repeats 3                                \
    --network_filename      networks/fsolv-${FFNAME}-network.json   \
    --scope_name_campaign   "${FFNAME_UNDERSCORED}"                 \
    --scopekey_output       scoped-keys/fsolv_${FFNAME}_key.dat     \
    --scope_name_project    "freesolv" > logs/submit/submit-freesolv-${FFNAME}.log
