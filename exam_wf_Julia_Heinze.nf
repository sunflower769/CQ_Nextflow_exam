nextflow.enable.dsl = 2

params.ref_accession = "M21012"
params.out = "$projectDir/output"
params.storeDir = "$projectDir/ref"
params.dataDir = "$projectDir/data"

process download_gb {
	storeDir params.storeDir
	input:
		val ref_accession
	output:
		path "${params.ref_accession}.fasta"
	"""
	wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${params.ref_accession}&rettype=fasta&retmode=text" -O ${params.ref_accession}.fasta
	"""	
}


process combine_fasta {
	publishDir params.out, mode:"copy", overwrite:true
	input:
		path ref_file
	output:
		path "ref_${params.ref_accession}_with_new_seqs.fasta"
	"""
	cat * >> ref_${params.ref_accession}_with_new_seqs.fasta
	"""
}

//potential_problem: line breaks in FASTA!!!

process mafft {
	publishDir params.out, mode:"copy", overwrite:true
	container "https://depot.galaxyproject.org/singularity/mafft:7.525--h031d066_1"
	input:
		path infile
	output:
		path "mafft_output_ref_${params.ref_accession}.fasta"
	"""
	mafft $infile > mafft_output_ref_${params.ref_accession}.fasta
	"""
}

process trimal {
	publishDir params.out, mode:"copy", overwrite:true
	container "https://depot.galaxyproject.org/singularity/trimal:1.5.0--h4ac6f70_1"
	input:
		path infile
	output:
		path "trimal_output_ref_${params.ref_accession}.fasta"
		path "trimal_report_ref_${params.ref_accession}.html"
	"""
	trimal -in $infile -out trimal_output_ref_${params.ref_accession}.fasta -fasta -automated1 -htmlout trimal_report_ref_${params.ref_accession}.html
	"""
}


workflow {
	refCh = download_gb(Channel.from(params.ref_accession)) 
	seqCh = Channel.fromPath("${params.dataDir}/*.fasta")
	allCh = refCh.concat(seqCh) | collect | combine_fasta | mafft | trimal
}