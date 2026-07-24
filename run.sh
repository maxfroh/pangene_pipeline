#!/bin/bash
ENV_NAME="pangene_pipeline"

if mamba env list | grep -q "^${ENV_NAME}\s"; then
    :
else
    echo "Error: Required environment \"${ENV_NAME}\" does not exist. Building now."
    mamba env create -f environment.yml -y

    # install downstream_analysis script
    KAKS_FILE="$CONDA_PREFIX/bin/add_ka_and_ks_to_collinearity.pl"
    wget https://raw.githubusercontent.com/wyp1125/MCScanX/refs/heads/master/downstream_analyses/add_ka_and_ks_to_collinearity.pl -o $KAKS_FILE -nv
    chmod +x $KAKS_FILE
fi

mamba run -n $ENV_NAME python pipeline.py "$@"