# Pipeline Configuration Files

To run this code, you must construct a [TOML](https://toml.io/en/) configuration file. A template is included.

There are a variety of required and non-required fields to fill out. 

TOML files have the following structure:
```toml
[section]
key1 = value1

[section.subsection]
key2 = value2
```
`section.key1` is equal to `value1`, while `section.subsection.key2` is equal to `value2`.

Single-line tables can look like this:
```toml
[section]
basic_table = {"key1": "val1", "key2": "val2"}
table_of_lists = {"key1": ["val1a", "val1b"], "key2": ["val2a"]}
```

## Configuration File Sections
|Section Name|Has Subsections?|Description|
|:-|:-:|:-|
|`input`|✗|This section of the configuration file contains all required general input information.|
|`output`|✗|This section specifies where results and output should be kept.|
|`cores`|✗|This section specifies how many CPU cores are available.|
|`parameters`|✗|This section contains all hyperparameters for the pipeline.|
|`pangene`|✓|This section contains at least one subsection containing information for a pangene reference. The subsection should be named appropriately (e.g., `pangene.PAN_1`, `pangene.PAN_2`, etc.)|
|`reference`|✓|This section contains at least one subsection containing the annotation and FASTA information for a single-reference genome. The subsection should be named appropriately (e.g., `reference.REF_1`, `reference.REF_2`, etc.)|
|`runs`|✓|This section contains at least one subsection indicating a single run of the pipeline. Each run contains a single group of samples to analyze. The subsection should be named appropriately (e.g., `runs.RUN_1`, `runs.RUN_2`, etc.)|

## Configuration File Keys
|Key|Location|Required?|Type|Default Value|Description|
|-|-|:-:|:-:|-|-|
|`base_dir`|`input`|✓|`string`|n/a|The base directory where data is kept and stored.|
|`results_dir`|`output`|✓|`string`|n/a|The directory where results will be stored.|
|`out_file`|`output`|✗|`string`|"out.txt"|The file where `stdout` will be redirected.|
|`err_file`|`output`|✗|`string`|"err.txt"|The file where `stderr` will be redirected.|
|`p`|`cores`|✗|`int`|1|The total number of processors available to the pipeline.|
|`auto_allocate_processors`|`cores`|✗|`boolean`|false|Whether to allow the pipeline to calculate the number of processors available. Will allow the tool to use all available processors!|
|`alpha`|`parameters`|✗|`float`|0.05|The value to use for p-value comparison ($p < \alpha$) when determining differential expression determination.|
|`l2FC_thresh`|`parameters`|✗|`float`|1|The value to use for absolute $log_2$-fold-change thresholding when determining differential expression.|
|`frag_length_mean`|`parameters`|✗|`float`|200|The mean read fragment length to provide to kallisto.|
|`frag_length_std`|`parameters`|✗|`float`|20|The standard deviation of read fragment length to provide to kallisto.|
|`grp_file`|`pangene.PANGENE_NAME`|✓|`str`|n/a|The melted `grp` file to reference for the pangene. May be replaced later.|
|`pangene_fastas_dir`|`pangene.PANGENE_NAME`|✓|`str`|n/a|The directory containing all FASTA files for the genomes the pangene is composed of. May be replaced later.|
|`redundancy_thresh`|`pangene.PANGENE_NAME`|✗|`float`|0.98|The sequence identity threshold for CD-HIT. See the [CD-HIT wiki](https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki#user-content-CDHITEST) for more information.|
|`annotation_file`|`reference.REFERENCE_NAME`|✓|`str`|n/a|The file to use for annotations for the reference `REFERENCE_NAME` (`(.gff\|.gtf\|.gff3)[.gz]` file expected).|
|`cds_fasta`|`reference.REFERENCE_NAME`|✓|`str`|n/a|The FASTA file to use for the reference `REFERENCE_NAME` (`.fa[.gz]` file expected).|
|`use`|`run.RUN_NAME`|✓|`table`|n/a|A table containing two keys: `pangene` and `reference`, where the value of `pangene` is a list of all pangenes to use and the value of `reference` is a list of all single-genome references to use. These names must be established in the configuration file (i.e., `pangene.PANGENE_NAME` or `reference.REFERENCE_NAME` must be defined).|
|`conditions`|`runs.RUN_NAME`|✓|`list` of `string`|n/a|The names of samples being tested. If there are not multiple replicates, this can be the name of the different samples themselves. Otherwise, it is expected they would follow a pattern where the sample name is a substring of the replicates' names. Example: `SampleName` and replicates `SampleName_R1.fq.gz`, `SampleName_R2.fq.gz`, etc.|
|`sample_dir`|`runs.RUN_NAME`|✓|`string`|n/a|The directory where sample files are stored.|
|`samples`|`runs.RUN_NAME`|✓|`list` of `string`|n/a|The sample files to use (`.fq[.gz]` files expected).|