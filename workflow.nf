#!/usr/bin/env nextflow

params.bamfile = '../../input/mt.bam'
params.bedfile = '../../input/regions.bed.gz'

process COUNTING {
	// TODO: Need to add bioconda to channel
	// conda "bioconda::samtools=1.14"
	
	input:
	path bamfile
	path bedfile

	output:
	stdout 

	script:
	"""
	samtools sort ${bamfile} -o sorted_${bamfile}
	samtools index sorted_${bamfile}
	samtools bedcov ${bedfile} sorted_${bamfile}
	"""
}

// TODO: Output needs to be JSON

process EXTRACT {
	
	input:
	path bedfile
	path bamfile

	output:
	path 

	script:
	"""

	"""
}

workflow{
	bamfile_ch = Channel.fromPath(params.bamfile)
	bedfile_ch = Channel.fromPath(params.bedfile)

	COUNTING(bamfile_ch,
		bedfile_ch)
}
