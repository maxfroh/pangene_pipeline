# pangene_pipeline


```mermaid
%%{init: {'theme': 'dark'}}%%
classDiagram

    class ParamManager {
        + p: int
        + alpha: float
        + l2FC_thresh: float
        + frag_length_mean: int
        + frag_length_std: int
    }
    class PangeneConstructor {
        + build_pangene(pm: ParamManager)
        + get_pangene_reference_manager()
    }
    class ReferenceManager {
        + reference: str
        + annotation_file : Path
        + cds_fasta: Path
        + tmp_dir: Path
        + kallisto_dir: Path
        + out_file: Path
        + err_file: Path
        + create_dirs()
    }
    class RunManager {
        + run: str
        + base_dir: Path
        + sample_dir: Path
        + results_dir: Path
        + out_base_file: Path | str
        + err_base_file: Path | str
        + param_manager: ParamManager
        + reference_managers: list[ReferenceManager]
        + run_deg()
        - create_dirs()
    }
    
    class input {
        + base_dir: string
    }
    class output {
        + results_dir: string
        + out_file_template: string
        + err_file_template: string
    }
    class parameters {
        + p: integer
        + auto_allocate_processors: boolean
        + alpha: float
        + l2FC_thresh: float
        + frag_length_mean: integer
        + frag_length_std: integer
    }
    class pangene {
        name
        + grp_file: string
        + pangene_fastas_dir: string
        + redundancy_thresh: float
    }
    class reference {
        name
        + annotation: string
        + cds: string
    }
    class run {
        name
        + use: table[string, array[string]]
        + conditions: array[string]
        + samples_dir: string
        + samples: array[string]
    }

    class config.toml
    config.toml "1" *-- "1" input
    config.toml "1" *-- "1" output
    config.toml "1" *-- "1" parameters
    config.toml "1" *-- "0..*" pangene
    config.toml "1" *-- "0..*" reference
    config.toml "1" *-- "1..*" run
    run "1" *-- "0..*" parameters : may contain
    RunManager "1" *-- "1..*" ReferenceManager : contains
    RunManager "1" *-- "1" ParamManager : contains
    RunManager ..|> config.toml : gets variables from
    ReferenceManager ..|> PangeneConstructor : can be created from

```