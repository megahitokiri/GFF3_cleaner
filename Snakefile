#--------------------------------------------------------------------------------
# MAIN CONFIGURATION PARAMETERS (INTERNAL CONFIG FILE)
#Cancel all jobs if fail: squeue -u $USER -h | awk '{print $1}' | xargs scancel
#--------------------------------------------------------------------------------

CHRS = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18'.split()
PROJECT = "HanIR"
REFERENCE = "MAIN_FASTAs/HanIRr1.0-20201123.genome.fasta"
NEW_NAMES_REFERENCE = "HanIRr1.0-20201123.genome.new_names.fasta"
GFF3_FILE = "MAIN_GFF3s/HanIRr1.0-20201123.gff3"

#--------------------------------------------------------------------------------
# TargetRule FINAL_GFF3
#--------------------------------------------------------------------------------

rule FINAL_GFF3:
	input:
		expand("{Main_Reference}",Main_Reference=REFERENCE),
		expand("{Project}/Ref/scaffold_{Chrs}.fasta",Project=PROJECT,Chrs = CHRS),
		expand("{Project}/EDTA_Files/scaffold_{Chrs}.fasta.mod.EDTA.TElib.fa",Project=PROJECT,Chrs = CHRS),
		expand("{Project}/EDTA_Files/scaffold_{Chrs}.fasta.mod.EDTA.intact.gff3",Project=PROJECT,Chrs = CHRS),
		expand("{Project}/EDTA_Files/scaffold_{Chrs}.masked.fasta",Project=PROJECT,Chrs = CHRS),
		expand("{Gff3_file}",Gff3_file=GFF3_FILE),		
#--------------------------------------------------------------------------------
# Init: Initializing files and folder
#--------------------------------------------------------------------------------
rule Init:
	input:
		reference={REFERENCE},
		Assembly_spliter="Assembly_Chr_splitter.R",
	output:
		expand("{Project}/Ref/{New_names_ref}",Project=PROJECT,New_names_ref=NEW_NAMES_REFERENCE),
	params:
		project=PROJECT,
	shell:
		"""
		snakemake --dag | dot -Tsvg > dag.svg
		mkdir {params.project}
		cd {params.project}
		
			mkdir logs
			mkdir EDTA_Files
			mkdir Ref
			mkdir Annotation_steps
			mkdir FINAL_ANNOTATION
			mkdir COGNATE
		cd ..
		cp  {input.reference} {output}

		ml samtools

		#Index FASTA file
		samtools faidx {output} 
		
		cd {params.project}
		
		ml unload samtools
		"""

#--------------------------------------------------------------------------------
# Chr_splitting: Split the Main Assembly in Chromosomes for easy handling
#--------------------------------------------------------------------------------

rule Chr_splitting:
	input:
		rules.Init.output,
	output:
		expand("{Project}/Ref/scaffold_{Chrs}.fasta",Project=PROJECT,Chrs=CHRS),
		
	params:
		project=PROJECT,
	shell:
		"""
		echo Assembly_Chr_splitter.R --args -f {input} 
		ml r/4.1.0

		echo Assembly split into Chromosomes
		R --vanilla < Assembly_Chr_splitter.R --args -f {input} &&
		mv *scaffold*.fasta {params.project}/Ref/

		ml unload r/4.1.0
		ml unload samtools

		"""

#--------------------------------------------------------------------------------
# EDTA_individual: Look for TE elements on individual Fasta Chr
#--------------------------------------------------------------------------------
rule EDTA_individual:
	input:
		rules.Chr_splitting.output,
	output:
		fa_file="{Project}/EDTA_Files/scaffold_{Chrs}.fasta.mod.EDTA.TElib.fa",
		gff3_file="{Project}/EDTA_Files/scaffold_{Chrs}.fasta.mod.EDTA.intact.gff3",
	params:
		project=PROJECT,

	shell:
		"""
		cp {params.project}/Ref/*scaffold*.* {params.project}/EDTA_Files

		cd {params.project}/EDTA_Files
		eval "$(conda shell.bash hook)"

		conda activate EDTA
			echo starting EDTA process on: scaffold_{wildcards.Chrs}.fasta
			EDTA.pl --genome scaffold_{wildcards.Chrs}.fasta
		conda deactivate
		"""
		
#--------------------------------------------------------------------------------
# Masked_FASTA: Create masked fasta for further analysis from EDTA results.
#--------------------------------------------------------------------------------

rule Masked_FASTA:
	input:
		EDTA_repeats_file=rules.EDTA_individual.output.gff3_file,
		reference=rules.Chr_splitting.output,
	output:
		masked_fasta_file="{Project}/EDTA_Files/scaffold_{Chrs}.masked.fasta",
	params:
		project=PROJECT,
	shell:
		"""
		ml bedtools
		cd {params.project}/EDTA_Files
		echo "Creating Masked Reference Genome for scaffold_{wildcards.Chrs}.masked.fasta"
		bedtools maskfasta -fi scaffold_{wildcards.Chrs}.fasta -bed scaffold_{wildcards.Chrs}.fasta.mod.EDTA.intact.gff3 -fo {params.project}_{wildcards.Chrs}.masked.fasta
		ml unload bedtools	
		"""

#--------------------------------------------------------------------------------
# STEP1_annotation: Transform Analysis .
#--------------------------------------------------------------------------------

rule STEP1_annotation:
	input:
		Masked_FASTA_file=rules.Masked_FASTA.output.masked_fasta_file,
		Gff3_file={GFF3_FILE},
	output:
		masked_fasta_file="{Project}/EDTA_Files/scaffold_{Chrs}.masked.fasta",
	params:
		project=PROJECT,
	shell:
		"""
		ml r/4.1.0
		cp {input.Gff3_file} {params.project}/Annotation_steps/{params.project}.gff3
		echo "Processing Step1 of annotation"
		R --vanilla < Step1_Annotation.R --args -a {params.project} -c {wildcards.Chrs}
		ml unload r/4.1.0
		"""
