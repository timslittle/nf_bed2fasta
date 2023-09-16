# nf_bed2fasta
Simple Nextflow pipeline to count read alignments, and extract sequences, from a bam file, corresponding to regions specified by an input bed file.

## Installation.

## Usage.

If your HPC system has a Nextflow profile then be sure to use this after the `-profile` flag too. You can specify multiple profiles as follows: `-profile test,docker`.

`nextflow run workflow.nf --bamfile <INPUT.bam> --bedfile <INPUT.bed> -profile docker`

It is recommended to use `-profile docker`, however you could omit this if you have `samtools` and its dependencies running locally.

After a successful run, the `results` folder will contain the counts file and sequence fasta file.

## Tests.

## Tips and tricks.

*   You can use `-w` to specify the path for the `work` directory. This directory can get very large with Nextflow and is usually best deleted after use, however you do need to keep it to `-resume` pipelines using the caching system. So, if you plan to run this pipeline multiple times, it may be worth specifying different `work` directories so that you can delete only the ones you definitely no longer require and free up a lot of space.
