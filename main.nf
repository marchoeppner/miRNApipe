#!/usr/bin/env nextflow

/*

||  ||
miRNApipe
||  ||

miRNApipe alignment and quantification pipeline
Concept: S. Juzenas
Implementation: M. Hoeppner & S. Juzenas

*/

try {
    if( ! nextflow.version.matches(">= $workflow.manifest.nextflowVersion") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $workflow.manifest.nextflowVersion required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please use a more recent version of Nextflow!\n" +
              "============================================================"
}

def helpMessage() {
  log.info"""
============================================
>> miRNApipe alignment and quantification <<
============================================

Usage:

A typical command to run this pipeline would be:

nextflow run main.nf --reads 'data/*_R{1,2}_001.fastq.gz' --assembly GRCh37 --kit Nextera

Mandatory arguments:

--reads                 A pattern to define which reads to use
--genome		The  name of the genome assembly to use
--kit			The library/adapter kit that was used (may alternatively specify the adapters as FASTA file, see below). 

Options:
--igenomes_base		Location of the iGenomes folder on your cluster (may also be specified in your cluster config file)
--adapters              A gzip compressed FASTA file with sequencing adapters (if no kit is specified!)
--email                 Provide an Email to which reports are send.
--run_name              A name for this run.
--skip_multiqc          Do not generate a QC report
""".stripIndent()
}


// Show help message
if (params.help){
        helpMessage()
        exit 0
}

log.info "=================================================="
log.info "miRNApipe alignment and quantification v${workflow.manifest.version}"
log.info "Nextflow Version:     $workflow.nextflow.version"
log.info "Command Line:         $workflow.commandLine"
log.info "Authors:              S. Juzenas & Marc Hoeppner"
log.info "================================================="
log.info "Starting at:          $workflow.start"

/* +++++++++++++++++++++++++++++++++
Specify all input settings and files
+++++++++++++++++++++++++++++++++ */

if (params.kit) {
	params.adapter = params.kits[params.kit].adapter
	params.adapter_2 = params.kits[params.kit].adapter_2
} else if (params.adapter) {

} else {
	// exit 1; "Don't know which adapter(s) to use for trimming (use --kit or --adapter/--adapter_2)"
}

params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false

if (params.star_index) {
	star_index = Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
} else if (params.fasta) {
	Channel.fromPath(params.fasta)
        .ifEmpty { exit 1, "Fasta file not found: ${params.fasta}" }
	.into { ch_fasta_for_star_index  }
} else {
    exit 1, "No reference genome specified!"
}

// Whether to send a notification upon workflow completion
params.email = false

if(params.email == false) {
	exit 1, "You must provide an Email address to which pipeline updates are send!"
}

def summary = [:]

OUTDIR=file(params.outdir)

run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

summary['runName'] = run_name
if (params.genome) {
	summary['genome'] = params.genome
}
if (params.fasta ) {
	summary['fasta'] = params.fasta
}
summary['gtf'] = params.gtf
summary['reads'] = params.reads
summary['Current home'] = "$HOME"
summary['Current user'] = "$USER"
summary['Current path'] = "$PWD"
summary['SessionID'] = workflow.sessionId

//
// PIPELINE STARTS HERE
//

Channel
	.fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
	.ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
	.into { raw_reads_fastqc; raw_reads_fastp }


// If no STAR index + resources are available

if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtf_makeSTARindex; gtf_star }
} else {
    exit 1, "No GTF annotation specified!"
}

if(!params.star_index && params.fasta){
   
	 process makeSTARindex {
		tag "$fasta"
		publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

		input:
		file fasta from ch_fasta_for_star_index
		file gtf from gtf_makeSTARindex

		output:
		file "star" into star_index

		script:
		def avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
		"""
		mkdir star
		STAR \\
		--runMode genomeGenerate \\
		--runThreadN ${task.cpus} \\
		--sjdbGTFfile $gtf \\
		--genomeDir star/ \\
		--genomeFastaFiles $fasta \\
		$avail_mem
		"""
	}
} else if (params.star_index) {
	star_index = Channel.fromPath(params.star_index)
} else {
	exit 1; "Neither a star index folder (--star_index) nor a genome fasta file (--fasta) provided, cannot proceed."
}
// Clean up reads

process runFastp {

	tag "${id}"
	publishDir "${OUTDIR}/${id}/fastp", mode: 'copy'

	input:
	set val(id),file(reads) from raw_reads_fastp

	output:
	set val(id),file("*_trimmed.fastq.gz") into trimmed_reads
	set file(html),file(json) into fastp_stats
	
	script:
	def options = ""
	
	if (params.adapter) {
		options += " -a ${params.adapter}"
	} 
	if (params.adapter && !params.singleEnd ) {
		options += " --adapter_sequence_r2 ${params.adapter_2}"
	}

        left = file(reads[0]).getName() + "_trimmed.fastq.gz"
	json = file(reads[0]).getBaseName() + ".fastp.json"
	html = file(reads[0]).getBaseName() + ".fastp.html"

	if (params.singleEnd) {
		"""
                fastp $options --in1 $reads --out1 $left -w ${task.cpus} -j $json -h $html --length_required 14 -p

		"""
	} else {
                right = file(reads[1]).getName() + "_trimmed.fastq.gz"

		"""
		fastp $options --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right -w ${task.cpus} -j $json -h $html --length_required 14 -p
		"""
	}

}

process runStar {

	tag "${id}"
        publishDir "${OUTDIR}/${id}/STAR", mode: 'copy'

	input:
	set val(id),file(reads) from trimmed_reads
        file index from star_index.collect()
        file gtf from gtf_star.collect()

	output:
	set val(id),file(bam),file(bai) into starBam
	set val(id),file(counts) into StarCounts
	set val(id),file(wig_plus),file(wig_minus) into StarWiggle

	script:
	bam = id + "Aligned.sortedByCoord.out.bam"
	bai = bam + ".bai"
	counts = id + "ReadsPerGene.out.tab"
	wig_plus = id + "Signal.Unique.str1.out.wig"
	wig_minus = id + "Signal.Unique.str2.out.wig"
		
	"""
		STAR --genomeDir ${params.star_index} --readFilesIn $reads  \
		--outSAMunmapped Within \
       		--readFilesCommand zcat \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode TranscriptomeSAM GeneCounts \
		--outReadsUnmapped Fastx \
		--alignEndsType EndToEnd \
		--outFilterMismatchNmax 1 \
		--outFilterMultimapScoreRange 0 \
		--outFilterScoreMinOverLread 0 \
		--outFilterMatchNminOverLread 0 \
		--outWigType wiggle \
       		--outFileNamePrefix $id \
		--sjdbGTFfile $gtf \
       		--limitBAMsortRAM 32212254720 \
       		--runThreadN ${task.cpus} --genomeLoad NoSharedMemory \
       		--outFilterMultimapNmax ${params.n_multimap} \
		--alignSJDBoverhangMin 1000 \
		--alignIntronMax 1 \
		--outWigStrand Stranded \
		--outWigNorm RPM \
		--clip3pAdapterSeq TGGAATTCTC --clip3pAdapterMMp 0.1 --outFilterMismatchNoverLmax 0.03

		samtools index $bam
	"""
}





