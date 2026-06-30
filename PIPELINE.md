# Pangene

* *TBD* - orthofinder, mcscanx, paml
    * Get melted `.grp` file
* Load `.grp` file via pandas to get all individual genes/groups and group by genome
* **&#8741; (8?)** for all CDS fastas:  `gunzip` genome fasta &rarr; `bgzip` to `.bgz` &rarr; `samtools faidx` on the bgz file
* create a dictionary (`defaultdict(list)`) to store fasta sequences for each OG 
* **&#8741; (8?)** for each genome: extract fasta sequences and add to appropriate dictionary entries
* write dictionary entries to fastas for each OG
* **&#8741; (max)** cd-hit-est on each OG fasta
* concatenate all files using `find _ > xargs cat`


# Kallisto
* *TBD?* clean reads with trimmomatic
* build a `kallisto index` with transcriptome reference
* for each sample, `kallisto quant`
* send kallisto results directly to `R`
    * read annotation file
    * tximport
    * DESeq2
* ?