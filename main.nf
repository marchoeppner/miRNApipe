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

OUTDIR=file(params.outdir)

if (params.kit) {
	params.adapter = params.kits[params.kit].adapter
	params.adapter_2 = params.kits[params.kit].adapter_2
} else if (params.adapter) {

} else {
	// exit 1; "Don't know which adapter(s) to use for trimming (use --kit or --adapter/--adapter_2)"
}

params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.ftp = params.genome ? params.genomes[ params.genome ].ftp ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false

if (!params.gtf) {
	exit 1, "Missing miRNA annotations in GTF/GFF format (--gtf)"
}

// Figure out where the genome reference comes from
// 1. We have a pre-existing STAR index
if (params.star_index) {
	star_index = Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
// Not a pre-existing index, but at least we have a fasta file
} else if (params.fasta) {
	println "Using user-provided FASTA file to build STAR index."
	Channel.fromPath(params.fasta)
        .ifEmpty { exit 1, "Fasta file not found: ${params.fasta}" }
	.into { ch_fasta_for_star_index  }

// If we only have a genome name, check if we at least have an FTP link to download it
} else if (params.ftp ) {

	println "No genome sequence or STAR index found, will try to download the correct assembly instead..."

	Channel.from(params.ftp)
		.set { link_to_fasta }

	process runFetchGenome {

                executor = "local"

		publishDir "${OUTDIR}/reference_genome", mode: 'copy'
			
		input:
		val(ftp_link) from link_to_fasta

		output:
		file(fasta) into ch_fasta_for_star_index

		script:

		fasta_gz = ftp_link.split("/")[-1]
		fasta = fasta_gz.replaceAll(".gz", "")

		"""
			wget $ftp_link

			gunzip $fasta_gz

		"""
	}

} else {
    exit 1, "No reference genome specified!"
}

// Whether to send a notification upon workflow completion
params.email = false

if(params.email == false) {
	exit 1, "You must provide an Email address to which pipeline updates are send!"
}

def summary = [:]

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
	.fromPath("$baseDir/assets/miRBase/v22/hairpin.fa.gz")
	.set { hairpin_for_index }

Channel
	.fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
	.ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
	.into { raw_reads_fastqc; raw_reads_fastp }


if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtf_makeSTARindex; gtf_star }
} else {
    exit 1, "No GTF annotation specified!"
}

if(!params.star_index && params.fasta || !params.star_index && params.ftp ){
   
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

process runHairpinStarIndex {

	publishDir "${OUTDIR}/hairpins", mode: 'copy'
	
	input:
	file(hairpins) from hairpin_for_index

	output:
	file(hairpin_index) into hairpin_index
	file(hairpin_fa) into hairpin_fasta

	script:
	hairpin_index = "star"
	hairpin_fa = hairpins.getBaseName()

	"""
	mkdir -p $hairpin_index
	
	gunzip -c $hairpins > $hairpin_fa

	STAR \
        	--runMode genomeGenerate \
                --runThreadN ${task.cpus} \
                --genomeDir $hairpin_index \
                --genomeFastaFiles $hairpin_fa
	"""
}

process runFastp {

	tag "${id}"
	publishDir "${OUTDIR}/${id}/fastp", mode: 'copy'

	input:
	set val(id),file(reads) from raw_reads_fastp

	output:
	set val(id),file("*_trimmed.fastq.gz") into trimmed_reads, trimmed_reads_hairpin
	file(json) into fastp_stats
	
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

process runStarHairpins {

	publishDir "${OUTDIR}/${id}/STAR_Hairpins", mode: 'copy'
	
	input:
	file(star_index) from hairpin_index.collect()
	set val(id),file(reads) from trimmed_reads_hairpin

	output:
	set file(id),file(bam) into bam_hairpin

	script:
	bam = id + "Aligned.sortedByCoord.out.bam"
        bai = bam + ".bai"
        counts = id + "ReadsPerGene.out.tab"
        wig_plus = id + "Signal.Unique.str1.out.wig"
        wig_minus = id + "Signal.Unique.str2.out.wig"
        log_file = id + "Log.final.out"

	"""
                STAR --genomeDir $star_index --readFilesIn $reads  \
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

process runMirTop {

	publishDir "${OUTDIR}/${id}/mirTop", mode: 'copy'

	input:
	set val(id),file(bam) from bam_hairpin
	file(hairpin_fa) from hairpin_fasta.collect()

	output:
	file(out) into mirtop_out
	
	script:

	out = id + "_out"

	"""
		mirtop gff -sps hsa --hairpin $hairpin_fa --gtf $params.gtf -o $out $bam
	"""
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
	file(log_file) into LogToMultiqc
	file(counts) into CountsToMultiqc

	script:
	bam = id + "Aligned.sortedByCoord.out.bam"
	bai = bam + ".bai"
	counts = id + "ReadsPerGene.out.tab"
	wig_plus = id + "Signal.Unique.str1.out.wig"
	wig_minus = id + "Signal.Unique.str2.out.wig"
	log_file = id + "Log.final.out"
		
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

process runMultiqc {

	publishDir "${OUTDIR}/MultiQC", mode: 'copy'

	input:
	file(json) from fastp_stats.collect()
	file(log) from LogToMultiqc.collect()
	file(counts) from CountsToMultiqc.collect()

	output:
	file("*.html") into multiqc_report

	script:

	"""
		cp $baseDir/assets/multiqc_config.yaml . 
		multiqc . 
	"""	
}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="

  def email_fields = [:]
  email_fields['version'] = workflow.manifest.version
  email_fields['session'] = workflow.sessionId
  email_fields['runName'] = run_name
  email_fields['success'] = workflow.success
  email_fields['dateStarted'] = workflow.start
  email_fields['dateComplete'] = workflow.complete
  email_fields['duration'] = workflow.duration
  email_fields['exitStatus'] = workflow.exitStatus
  email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
  email_fields['errorReport'] = (workflow.errorReport ?: 'None')
  email_fields['commandLine'] = workflow.commandLine
  email_fields['projectDir'] = workflow.projectDir
  email_fields['script_file'] = workflow.scriptFile
  email_fields['launchDir'] = workflow.launchDir
  email_fields['user'] = workflow.userName
  email_fields['Pipeline script hash ID'] = workflow.scriptId
  email_fields['genome'] = REF
  email_fields['manifest'] = workflow.manifest
  email_fields['summary'] = summary

  email_info = ""
  for (s in email_fields) {
	email_info += "\n${s.key}: ${s.value}"
  }

  def output_d = new File( "${params.outdir}/pipeline_info/" )
  if( !output_d.exists() ) {
      output_d.mkdirs()
  }

  def output_tf = new File( output_d, "pipeline_report.txt" )
  output_tf.withWriter { w -> w << email_info }	

 // make txt template
  def engine = new groovy.text.GStringTemplateEngine()

  def tf = new File("$baseDir/assets/email_template.txt")
  def txt_template = engine.createTemplate(tf).make(email_fields)
  def email_txt = txt_template.toString()

  // make email template
  def hf = new File("$baseDir/assets/email_template.html")
  def html_template = engine.createTemplate(hf).make(email_fields)
  def email_html = html_template.toString()
  
  def subject = "Diagnostic exome analysis finished ($run_name)."

  if (params.email) {

  	def mqc_report = null
  	try {
        	if (workflow.success && !params.skip_multiqc) {
            		mqc_report = multiqc_report.getVal()
            		if (mqc_report.getClass() == ArrayList){
                		log.warn "[IKMB miRNA-pipeline] Found multiple reports from process 'multiqc', will use only one"
                		mqc_report = mqc_report[0]
                	}
        	}
    	} catch (all) {
        	log.warn "[IKMB ExoSeq] Could not attach MultiQC report to summary email"
  	}

	def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
	def sf = new File("$baseDir/assets/sendmail_template.txt")	
    	def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    	def sendmail_html = sendmail_template.toString()

	try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
        }

  }

}

