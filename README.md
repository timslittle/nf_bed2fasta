# nf_bed2fasta
Simple Nextflow pipeline to count read alignments, and extract sequences, from a bam file, corresponding to regions specified by an input bed file.

## Installation.

Clone this repo to your device using `git clone https://github.com/timslittle/nf_bed2fasta`. Alternatively, just download the `workflow.nf` and `nextflow.config` files.

## Dependencies

In order to work, this pipeline requires Nextflow to be installed. Please visit https://www.nextflow.io/docs/latest/getstarted.html#installation for the latest installation instructions.

At a minimum, the pipeline also requires one of the following:
*   Locally installed `samtools` and `R`.
*   Docker.
*   Singularity.
*   (Ana)conda.

## Usage.

`nextflow run workflow.nf --bamfile <INPUT.bam> --bedfile <INPUT.bed> -profile docker`

If your HPC system has a Nextflow profile then be sure to use this after the `-profile` flag too. You can specify multiple profiles as follows: `-profile test,docker`.

It is recommended to use a `-profile` such as `docker`, `singularity` or `conda`. However you could omit specifying a profile completely if you have `samtools`, `R`, and their own dependencies running locally.

After a successful run, the `results` folder will contain the counts file inside `COUNTS` and sequence fasta file inside `EXTRACT'.

## Tests.

_Under development_. Use `-profile test` to run the pipeline using a minimal dataset, and test for the expected output by searching for sequences we know should be returned, and md5 checksums.

## Tips and tricks.

*   You can use `-w` to specify the path for the `work` directory. This directory can get very large with Nextflow and is usually best deleted after use, however you do need to keep it to `-resume` pipelines using the caching system. So, if you plan to run this pipeline multiple times, it may be worth specifying different `work` directories so that you can delete only the ones you definitely no longer require and free up a lot of space.
