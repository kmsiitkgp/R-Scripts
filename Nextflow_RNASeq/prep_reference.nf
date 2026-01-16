process PREP_REFERENCE {

	input:
	
	
	script:
	
	"""
	
	# 1. Create STAR index
	mkdir -p star_index_dir
	STAR --runMode genomeGenerate \\
        --runThreadN "${task.cpus}" \\
        --genomeDir star_index_dir \\
        --genomeFastaFiles "${params.ref_fasta}" \\
        --sjdbGTFfile "${params.ref_gtf}" \\
        --sjdbOverhang 100 \\
        --genomeSAindexNbases 14 \\
		1>> "STAR.error.log" 2>&1 \\
		&& echo "✅ STAR index generation completed successfully." \\
		|| { echo "❌ STAR index generation failed. Check STAR.error.log"; exit 1; }
		
	# 2. Create BED of genome for RSEQC
	gtfToGenePred "${params.ref_gtf}" tmp_file
	genePredToBed tmp_file rseqc_bed \\
		&& echo "✅ RSEQC BED generation completed successfully." \\
		|| { echo "❌ RSEQC BED generation failed. Check RSEQC.error.log"; exit 1; }

	# 3. Create SALMON index
	# Create a decoy list from genome
	grep "^>" "${params.ref_fasta}" | cut -d " " -f1 | sed 's/>//' > salmon_decoy
	
	# Prepare transcriptome
	gffread "${params.ref_gtf}" \\
		-g "${params.ref_fasta}" \\
		-w salmon_transcriptome_fasta	

	# Prepare gentrome (transcripts+genome)
	cat salmon_transcriptome_fasta "${params.ref_fasta}" > salmon_gentrome_fasta
	
	# Build index	
    salmon index \\
        --transcripts salmon_gentrome_fasta \\
		--decoys salmon_decoy \\
		--threads "${task.cpus}" \\
		--kmerLen 31 \\
		--index salmon_index_dir \\
		1>> "SALMON.error.log" 2>&1 \\
		&& echo "✅ SALMON index generation completed successfully." \\
		|| { echo "❌ SALMON index generation failed. Check SALMON.error.log"; exit 1; }
	
	
	
	
	
	"""
	
	output:
	path("star_index_dir"), 	emit: star_index_dir
	path("rseqc_bed"), 			emit: rseqc_bed
	path("salmon_index_dir"), 	emit: salmon_index_dir
	
}