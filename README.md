# pangene_pipeline

## Use
This tool requires an environment manager. It is recommended you use [mamba](https://mamba.readthedocs.io/en/latest/#), though [conda](https://docs.conda.io/en/latest/) also works.


Once you have created your configuration file and installed the necessary dependencies, the pipeline can be run using the following command:
```bash
$ python pipeline.py --config config.toml
```

**Note:** in the future, some dependencies may be installed/found via a separate process and errors handled during installation/runtime. Currently, it will be assumed all external tools are on the user's `PATH`.


## Pipeline
```mermaid
flowchart LR
peptide_fastas@{ shape: docs }
gffs@{ shape: docs }
cds_fastas@{ shape: docs }
blastp_results@{ shape: docs }
mcscanx_results@{ shape: docs }
orthogroups@{ shape: doc }
downstream_results@{ shape: docs }
orthologs@{ shape: doc }
final@{ shape: doc }
blastp@{ shape: rounded }
mcscanx@{ shape: rounded }
orthofinder@{ shape: rounded }
mcscanx_downstream@{ shape: rounded }
subdivide_orthogroups@{ shape: rounded }
cd-hit@{ shape: rounded }

peptide_fastas-->blastp
blastp-->blastp_results
blastp_results-->mcscanx
gffs-->mcscanx
mcscanx-->mcscanx_results
mcscanx_results-->mcscanx_downstream
mcscanx_downstream-->downstream_results
cds_fastas-->orthofinder
orthofinder-->orthogroups
orthogroups-->subdivide_orthogroups
downstream_results-->subdivide_orthogroups
subdivide_orthogroups-->orthologs
orthologs-->cd-hit
cd-hit-->final
```
