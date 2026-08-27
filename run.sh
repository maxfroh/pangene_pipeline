#!/usr/bin/env bash
PANGENE_PIPELINE_ENV_NAME="pangene_pipeline"

if mamba env list | grep -q "^${PANGENE_PIPELINE_ENV_NAME}\s"; then
    :
else
    echo "Error: Required environment \"${PANGENE_PIPELINE_ENV_NAME}\" does not exist. Building now."
    mamba env create -f environment.yml -y

    ENV_LOC=$(mamba env list | grep -E "^\s*${PANGENE_PIPELINE_ENV_NAME}\s+" | awk -v env_name="$PANGENE_PIPELINE_ENV_NAME" '{ if ($2 ~ env_name) print $2; else print $3 }')

    # install downstream_analysis script
    KAKS_FILE="${ENV_LOC}/lib/add_ka_and_ks_to_collinearity.pl"
    wget -nv -O "$KAKS_FILE" https://raw.githubusercontent.com/wyp1125/MCScanX/refs/heads/master/downstream_analyses/add_ka_and_ks_to_collinearity.pl
    
    KAKS_EXE="${ENV_LOC}/bin/add_ka_and_ks_to_collinearity"
    echo "" | awk -v kaks_file=$KAKS_FILE '{ print "perl \""kaks_file"\" \"$@\"" }' > $KAKS_EXE
    chmod +x $KAKS_EXE
fi

mamba run -n "$PANGENE_PIPELINE_ENV_NAME" python pipeline.py "$@"
