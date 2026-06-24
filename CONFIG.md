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

## Configuration File Sections
|Section Name|Has Subsections?|Description|
|:-|:-:|:-|
|`input`|✗|This section of the configuration file will contain all required general input information.|
|`output`|✗|This section will specify where results and output should be kept.|
|`parameters`|✗|This section contains all hyperparameters for the pipeline.|
|`runs`|✓|This section will contain at least one subsection indicating a single run of the pipeline. Each run contains a single group of samples to analyze. The subsection should be named appropriately (e.g., `runs.RUN_1`, `runs.RUN_2`, etc.)|

## Configuration File Keys
|Key|Location|Required?|Type|Default Value|Description|
|-|-|:-:|:-:|-|-|
|`base_dir`|`input`|✓|`string`|n/a|The base directory where data is kept and stored.|
|`sample_dir`|`input`|✓|`string`|n/a|The directory containing samples.| 
|`results_dir`|`output`|✓|`string`|n/a|The directory where results will be stored.|
|`out_file`|`output`|✗|`string`|"out.txt"|The file where `stdout` will be redirected.|
|`err_file`|`output`|✗|`string`|"err.txt"|The file where `stderr` will be redirected.|
|`p`|`parameters`|✗|`int`|1|The number of processors available to the pipeline.|
|`auto_allocate_processors`|`parameters`|✗|`boolean`|false|Whether to allow the pipeline to calculate the number of processors available. Will allow the tool to use all available processors!|
|`alpha`|`parameters`|✗|`float`|0.05|The value to use for p-value comparison ($p < \alpha$) when determining differential expression determination.|
|`l2FC_thresh`|`parameters`|✗|`float`|1|The value to use for $log_2$-fold-change thresholding when determining differential expression.|
|`frag_length_mean`|`parameters`|✗|`float`|200|The mean read fragment length to provide to kallisto.|
|`frag_length_std`|`parameters`|✗|`float`|20|The standard deviation of read fragment length to provide to kallisto.|
|`annotations`|`runs.RUN_NAME`|?|`table`|n/a|A dictionary where the key is the reference name and the value is its annotation file (`.gff` or `.gtf`). Example: `{Ref1 = "ref1.gff3", Ref2 = "ref2.gff3"}`.|
|`references`|`runs.RUN_NAME`|✓|`table`|n/a|A dictionary where the key is the reference name and the value is its coding sequence file (`.fa[.gz]`). Example: `{Ref1 = "ref1_cds.fa.gz", Ref2 = "ref2_cds.fa.gz"}`.|
|`conditions`|`runs.RUN_NAME`|✓|`list` of `string`|n/a|The names of samples being tested. If there are not multiple replicates, this can be the name of the different samples themselves. Otherwise, it is expected they would follow a pattern where the sample name is a substring of the replicates' names. Example: `SampleName` and replicates `SampleName_R1.fq.gz`, `SampleName_R2.fq.gz`, etc.|
|`samples`|`runs.RUN_NAME`|✓|`list` of `string`|n/a|The sample files to use (`.fq[.gz]` files expected).|