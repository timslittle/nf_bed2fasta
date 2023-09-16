#!/usr/bin/env nextflow

params.bamfile = '../../input/mt.bam'
params.bedfile = '../../input/regions.bed.gz'

process COUNTING {
	
	input:
	path bedfile
	path bamfile

	output:
	path 

	script:
	"""

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
	COUNTING()
}
