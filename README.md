# pangene_pipeline

## Use
Once you have created your configuration file and installed the necessary dependencies, the pipeline can be run using the following command:
```bash
$ python pipeline.py --config config.toml
```

**Note:** in the future, some dependencies may be installed/found via a separate process and errors handled during installation/runtime. Currently, it will be assumed all external tools are on the user's `PATH`.


## Pipeline
```mermaid
flowchart LR
    pipeline["`This is the pipeline `"]
    node[s]
    pipeline-->node
```
