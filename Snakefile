#--------------------------------------------------------------------------------
# MAIN CONFIGURATION PARAMETERS (INTERNAL CONFIG FILE)
#Cancel all jobs if fail: squeue -u $USER -h | awk '{print $1}' | xargs scancel
#--------------------------------------------------------------------------------

CHRS = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18'.split()
PROJECT = "HanIR"
REFERENCE = "MAIN_FASTAs/HanIRr1.0-20201123.genome.fasta"
NEW_NAMES_REFERENCE = "HanIRr1.0-20201123.genome.new_names.fasta"

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
#		expand("Post_Maker_Files/MAKER_ORF_Filtered_{Project}.scaffold_{Chrs}.AED_{AED_filter}.gff3",Project=PROJECT,Chrs = CHRS,AED_filter=AED_FILTER),
		
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
	shell:
		"""
		ml bedtools
		cd EDTA_Files
		echo "Creating Masked Reference Genome for scaffold_{wildcards.Chrs}.masked.fasta"
		bedtools maskfasta -fi scaffold_{wildcards.Chrs}.fasta -bed scaffold_{wildcards.Chrs}.fasta.mod.EDTA.intact.gff3 -fo scaffold_{wildcards.Chrs}.masked.fasta
		ml unload bedtools	
		"""

