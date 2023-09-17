#!/usr/bin/env nextflow

params.bamfile = '../../input/mt.bam'
params.bedfile = '../../input/regions.bed.gz'

process COUNTING {
	// TODO: Need to add bioconda to channel
	// conda "bioconda::samtools=1.14"

	publishDir "results/${task.process}", mode: 'copy'

	// TODO: Ignoring this container command at the moment.
	// container 'quay.io/biocontainers/samtools:1.17--hd87286a_1'
	
	input:
	path bamfile
	path bedfile

	output:
	path "counts.json"

	script:
	"""
	samtools sort ${bamfile} -o sorted_${bamfile}
	samtools index sorted_${bamfile}
	samtools bedcov ${bedfile} sorted_${bamfile} | \\
		awk 'BEGIN{ print "\\{" } {print "\\"region\\":" "\\{\\"locus\\":""\\""\$1 ":" \$2 "-" \$3"\\"" ",\\"counts\\":" "\\""\$4"\\"\\}"} END{print "\\}"}' \\
		> counts.json
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
	bamfile_ch = Channel.fromPath(params.bamfile, checkIfExists: true)
	bedfile_ch = Channel.fromPath(params.bedfile, checkIfExists: true)

	COUNTING(bamfile_ch,
		bedfile_ch)
}
