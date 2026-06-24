_DEFAULT_PARAMS = {
    "p": 1,                 # num processers allocated
    # automatically determine the maximum number of processors available
    "auto_allocate_processors": False,
    "alpha": 0.05,          # for p-value filtering
    "l2FC_thresh": 1,        # for l2FC filtering
    "frag_length_mean": 200,  # for kallisto quant
    "frag_length_std": 20,   # for kallisto quant
}
"""
`DEFAULTS` contains default parameters to use for calculations, in case some are not
specified in the used configuration file. 
"""
