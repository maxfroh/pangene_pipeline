_DEFAULT_PARAMS = {
    "p": 1,  # the total number of cores allocated to the pipeline
    "auto_allocate_processors": False,  # automatically determine the maximum number of processors available
    "alpha": 0.05,  # for p-value filtering
    "l2FC_thresh": 1,  # for l2FC filtering
    "frag_length_mean": 200,  # for kallisto quant
    "frag_length_std": 20,  # for kallisto quant
    "redundancy_thresh": 0.98,  # by how much to reduce redunandant genes in the pangene
}
"""
`_DEFAULT_PARAMS` contains default parameters to use for calculations, in case some are not
specified in the used configuration file. 
"""

_ALL_STEPS = []

_OPTIMAL_CORES = {
    "internal": 8,
    "OrthoFinder_t": -1,  # -1 = max available
    "OrthoFinder_a": 8,
    "makeblastdb": 1,
    "blastn": 4,
    "MCScanX": 1,
    "CD-HIT": 8,
    "Trimmomatic": 8,
    "kallisto": 8,
}
