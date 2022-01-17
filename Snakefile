#--------------------------------------------------------------------------------
# MAIN CONFIGURATION PARAMETERS (INTERNAL CONFIG FILE)
#Cancel all jobs if fail: squeue -u $USER -h | awk '{print $1}' | xargs scancel
#--------------------------------------------------------------------------------

CHRS = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18'.split()
PROJECT = "HanHA300"
REFERENCE = "MAIN_FASTAs/HanHA300r0.9-20211022.genome.fasta"
NEW_NAMES_REFERENCE = "HanHA300r0.9-20211022.genome.new_names.fasta"
GFF3_FILE = "MAIN_GFF3s/HanHA300r0.9-20211022.gff3"

#--------------------------------------------------------------------------------
# TargetRule FINAL_GFF3
#--------------------------------------------------------------------------------

rule FINAL_GFF3:
	input:
		expand("{Main_Reference}",Main_Reference=REFERENCE),
		expand("{Project}/Ref/scaffold_{Chrs}.fasta",Project=PROJECT,Chrs = CHRS),
		expand("{Project}/EDTA_Files/scaffold_{Chrs}.fasta.mod.EDTA.intact.gff3",Project=PROJECT,Chrs = CHRS),
		expand("{Project}/EDTA_Files/{Project}_chr{Chrs}.masked.fasta",Project=PROJECT,Chrs = CHRS),
		expand("{Gff3_file}",Gff3_file=GFF3_FILE),
		expand("{Project}/Annotation_steps/{Project}_step1_chr{Chrs}.gff3",Project=PROJECT,Chrs = CHRS),
		expand("{Project}/Annotation_steps/{Project}_step2_chr{Chrs}.gff3",Project=PROJECT,Chrs = CHRS),
		expand("{Project}/FINAL_ANNOTATION/FINAL_{Project}_v1_1.sorted.gff3",Project=PROJECT),
		expand("{Project}/Summary_data/{Project}.protein.fasta",Project=PROJECT),
		expand("{Project}/Summary_data/JML.{Project}_GFF3_summary.txt",Project=PROJECT),
		expand("{Project}/Summary_data/Original.{Project}_GFF3_summary.txt",Project=PROJECT),
		expand("{Project}/Summary_data/{Project}_busco/short_summary.specific.eudicots_odb10.{Project}_busco.txt",Project=PROJECT),
	params:
		project=PROJECT,
	shell:
		"""
		echo Pipeline Finished correctly for {params.project}..... CONGRATS!!!
		echo Pipeline Finished correctly for {params.project}..... CONGRATS!!! > {params.project}.Pipeline_complete.txt
		mv *.err {params.project}/logs
		mv *.out {params.project}/logs

		date "+DATE: %D%nTIME: %T" >> {params.project}.Pipeline_complete.txt
		"""	
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
			mkdir Summary_data
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
		gff3_file="{Project}/EDTA_Files/scaffold_{Chrs}.fasta.mod.EDTA.intact.gff3",
	params:
		project=PROJECT,

	shell:
		"""
		cp {params.project}/Ref/*scaffold*.* {params.project}/EDTA_Files
		cp -v Name_checker_pre.sh {params.project}/EDTA_Files
		cp -v Name_checker_post.sh {params.project}/EDTA_Files
		
		cd {params.project}/EDTA_Files
		
		echo "Name_cheker pre and post correct EDTA bigger tha 15 characters name error on FASTA."
		bash Name_checker_pre.sh scaffold_{wildcards.Chrs}.fasta {params.project}
		
		eval "$(conda shell.bash hook)"
		conda activate EDTA
			echo starting EDTA process on: scaffold_{wildcards.Chrs}.fasta
			EDTA.pl --genome scaffold_{wildcards.Chrs}.fasta
		conda deactivate
		
		bash Name_checker_post.sh scaffold_{wildcards.Chrs}.fasta {params.project} scaffold_{wildcards.Chrs}.fasta.mod.EDTA.intact.gff3
		"""
		
#--------------------------------------------------------------------------------
# Masked_FASTA: Create masked fasta for further analysis from EDTA results.
#--------------------------------------------------------------------------------

rule Masked_FASTA:
	input:
		EDTA_repeats_file=rules.EDTA_individual.output.gff3_file,
		reference=rules.Chr_splitting.output,
	output:
		masked_fasta_file="{Project}/EDTA_Files/{Project}_chr{Chrs}.masked.fasta",
	params:
		project=PROJECT,
	shell:
		"""
		ml bedtools
		cd {params.project}/EDTA_Files
		echo "Creating Masked Reference Genome for scaffold_{wildcards.Chrs}.masked.fasta"
		bedtools maskfasta -fi scaffold_{wildcards.Chrs}.fasta -bed scaffold_{wildcards.Chrs}.fasta.mod.EDTA.intact.gff3 -fo {params.project}_chr{wildcards.Chrs}.masked.fasta
		ml unload bedtools	
		"""

#--------------------------------------------------------------------------------
# STEP1_annotation: Transform Analysis into parsed GFF3 files .
#--------------------------------------------------------------------------------

rule STEP1_annotation:
	input:
		Masked_FASTA_file=rules.Masked_FASTA.output.masked_fasta_file,
		Gff3_file={GFF3_FILE},
	output:
		gff3_step1="{Project}/Annotation_steps/{Project}_step1_chr{Chrs}.gff3",
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
#------------------------------------------------------------------------------------
# STEP2_annotation: Uses EDTA maked fasta file and ORF analysis to determine viable genes
#------------------------------------------------------------------------------------

rule STEP2_annotation:
	input:
		GFF3_File=rules.STEP1_annotation.output,
		Ref_File=rules.Chr_splitting.output,
		Masked_FASTA_File=rules.Masked_FASTA.output,
		Original_Gff3_file={GFF3_FILE},
	output:
		gff3_step2="{Project}/Annotation_steps/{Project}_step2_chr{Chrs}.gff3",
	params:
		project=PROJECT,
	shell:
		"""
		BASEDIR=$PWD 
		pwd
		cp -v $BASEDIR/Step2_Filtering.R {params.project}/Annotation_steps/
		cd {params.project}/Annotation_steps/
		ml r/4.1.0
		R --vanilla < Step2_Filtering.R --args -g $BASEDIR/{params.project}/Annotation_steps/{params.project}_step1_chr{wildcards.Chrs}.gff3 -a $BASEDIR/{params.project}/Ref/scaffold_{wildcards.Chrs}.fasta -m $BASEDIR/{params.project}/EDTA_Files/{params.project}_chr{wildcards.Chrs}.masked.fasta -o {params.project}_step2_chr{wildcards.Chrs} -s $BASEDIR/{input.Original_Gff3_file}
		ml unload r/4.1.0
		cd ..
		"""
		
#------------------------------------------------------------------------------------
# Chr_merge: Fuse all gff3 individual chromosomes into complete assembly again
#------------------------------------------------------------------------------------

rule Chr_merge:
	input:
		expand("{Project}/Annotation_steps/{Project}_step2_chr{Chrs}.gff3",Project=PROJECT,Chrs = CHRS),
	output:
		"{Project}/FINAL_ANNOTATION/FINAL_{Project}_v1_1.sorted.gff3",
	params:
		project=PROJECT,
		Chrs=CHRS,
	shell:
		"""
		cat {params.project}/Annotation_steps/{params.project}_step2_chr1.gff3 > {params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1.gff3
		
		for i in {{2..18}}
			do
			tail -n +4 {params.project}/Annotation_steps/{params.project}_step2_chr$i.gff3 >> {params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1.gff3
			done
		
		/home/jmlazaro/github/gff3sort/gff3sort.pl --chr_order original {params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1.gff3 > {params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1.sorted.gff3
	
		echo awk '{{gsub("character\\\(0\\\)", '0');print}}' {params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1.sorted.gff3 
		awk '{{gsub("character\\\(0\\\)", "0");print}}' {params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1.sorted.gff3 > {params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1_1.sorted.gff3
		"""		

#------------------------------------------------------------------------------------
# Summary_statistics: Get summary of the new gff3 file
#------------------------------------------------------------------------------------

rule Summary_statistics:
	input:
		GFF3_file=rules.Chr_merge.output,
		Ref_file=rules.Init.output,
	output:
		Protein_FASTA="{Project}/Summary_data/{Project}.protein.fasta",
	params:
		project=PROJECT,
	shell:
		"""
		BASEDIR=$PWD
		ml StdEnv/2020  gcc/9.3.0
		ml nixpkgs/16.09  gcc/5.4.0
		ml nixpkgs/16.09  gcc/7.3.0
		ml transdecoder/5.5.0
		
		gff3_file_to_proteins.pl --gff3 {input.GFF3_file} --fasta $BASEDIR/{input.Ref_file} --seqType prot > $BASEDIR/{params.project}/Summary_data/{params.project}.protein.fasta
		ml unload transdecoder/5.5.0
		
		echo Protein File Generated correctly CORRECTLY ....
		ml unload perl
		"""
#------------------------------------------------------------------------------------
# GFF3_statistics: Calculate Statistics for Original and filtered GFF3 files
#------------------------------------------------------------------------------------

rule GFF3_statistics:
	input:
		GFF3_file=rules.Chr_merge.output,
		Original_Gff3_file={GFF3_FILE},
	output:
		New_GFF3_summary="{Project}/Summary_data/JML.{Project}_GFF3_summary.txt",
		Original_GFF3_summary="{Project}/Summary_data/Original.{Project}_GFF3_summary.txt",
	params:
		project=PROJECT,
	shell:
		"""
		BASEDIR=$PWD
		cp -v GFF3_Summary_Statistics.R {params.project}/Summary_data/
		cd {params.project}/Summary_data/
		
		ml r/4.1.0
		echo Calculating Statistics on New processed file
		R --vanilla < GFF3_Summary_Statistics.R --args --gff $BASEDIR/{input.GFF3_file} -o JML.{params.project}
		
		echo Calculating Statistics on Original file
		R --vanilla < GFF3_Summary_Statistics.R --args --gff $BASEDIR/{input.Original_Gff3_file} -o Original.{params.project}
		cd ../..
		"""
		
#------------------------------------------------------------------------------------
# BUSCO: Evaluate the Cognate results into BUSCO protein mode
#------------------------------------------------------------------------------------

rule BUSCO:
	input:
		Protein_fasta=rules.Summary_statistics.output.Protein_FASTA,
	output:
		"{Project}/Summary_data/{Project}_busco/short_summary.specific.eudicots_odb10.{Project}_busco.txt",
	params:
		project=PROJECT,
	shell:
		"""
		BASEDIR=$PWD
		ml StdEnv/2020
		ml gcc/9.3.0
		ml openmpi/4.0.3
		ml busco/5.2.2
		mkdir {params.project}/Summary_data/busco_downloads
		mkdir {params.project}/Summary_data/busco_downloads/lineages
		
		cp /home/jmlazaro/BUSCO/eudicots_odb10.tar.gz {params.project}/Summary_data/busco_downloads/lineages
		cd {params.project}/Summary_data/busco_downloads/lineages
		tar -xvf eudicots_odb10.tar.gz
		cd $BASEDIR/{params.project}/Summary_data/
		
		busco -f -c 4 -m protein -i $BASEDIR/{params.project}/Summary_data/{params.project}.protein.fasta -o {params.project}_busco -l eudicots_odb10 --offline --download_path $BASEDIR/{params.project}/Summary_data/busco_downloads
		echo done BUSCO analysis
		"""
