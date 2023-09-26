#!/usr/bin/env nextflow

params.bamfile = '../input/mt.bam'
params.bedfile = '../input/regions.bed.gz'

process COUNTING {
	conda "bioconda::samtools=1.14"
	container 'quay.io/biocontainers/samtools:1.17--hd87286a_1'

	publishDir "results/${task.process}", mode: 'copy'
	
	input:
	path bamfile
	path bedfile

	output:
	path "counts.json"

	script:
	"""
	# Sort and index the .bam file
	samtools sort ${bamfile} -o sorted_${bamfile}
	samtools index sorted_${bamfile}

	# Use view to count the alignments within the regions specified by the bed file
	samtools view -hcM -L ${bedfile} sorted_${bamfile} > counts.txt

	# Echo the bed file and the counts file so awk will read them in as a single line.
	# 	Use awk to output as .json by specifying all the quotations, squiggly brackets, and colons. 
	#	The former two need to be escaped with backslash.
	echo \$(cat ${bedfile}) \$(cat counts.txt) | \\
		awk 'BEGIN{ print "\\{" } {print "\\"region\\":" "\\{\\"locus\\":""\\""\$1 ":" \$2 "-" \$3"\\"" ",\\"counts\\":" "\\""\$4"\\"\\}"} END{print "\\}"}' \\
		> counts.json
	"""
}

process BAM_2_SAM {

	conda "bioconda::samtools=1.14"
	container 'quay.io/biocontainers/samtools:1.17--hd87286a_1'

	input:
	path bamfile

	output:
	path "*.sam"

	script:
	"""
	#Get the prefix to save our files
	prefix=\$(echo ${bamfile} | sed -E "s/.bam\$//")
	samtools view -h -o \$prefix.sam ${bamfile}
	"""
}

process EXTRACT {

	publishDir "results/${task.process}", mode: 'copy'
	
	input:
	path samfile
	path bedfile

	output:
	path '*fasta'

	script:
	"""
	#!/usr/bin/env Rscript

		lsam <- readLines('mt.sam')
		dfsam <- lsam[!grepl(lsam, pattern = '^@', perl = TRUE)]
		sam <- read.table(text = dfsam, 
						fill = TRUE, 
						row.names = NULL,
						header = FALSE)

		# Subset sam by the bed file regions. 
		bed <- read.table('${bedfile}', 
						header = FALSE, 
						col.names = c('chr', 'start', 'end'))
		sam_chr <- sam[sam\$V3 == bed\$chr,]
		# Need to include those that start outside the region but overlap with it.
		index <- as.integer(sam_chr\$V4) >= as.integer(bed\$start) & 
		as.integer(sam_chr\$V4) <= as.integer(bed\$end) |
		as.integer(sam_chr\$V4) < as.integer(bed\$start) &
		as.integer(sam_chr\$V4) + nchar(sam_chr\$V10) >= as.integer(bed\$start)
		sam_reg <- sam_chr[index,]

		# Use pastes and cat to format the sequence names and sequences as a fasta file,
		#   as well as saving the output.
		cat(
		paste(
			apply(
			sam_reg, 
			1, 
			function(x) {
				paste( paste0('>', x[1]), x[10], sep = '\n')
			}),
			collapse = '\n'),
		file = 'mt.fasta')
	"""
}

workflow{
	bamfile_ch = Channel.fromPath(params.bamfile, checkIfExists: true)
	bedfile_ch = Channel.fromPath(params.bedfile, checkIfExists: true)

	COUNTING(bamfile_ch,
		bedfile_ch)
	
	samfile_ch = BAM_2_SAM(bamfile_ch)

	EXTRACT(samfile_ch,
	bedfile_ch)
}
