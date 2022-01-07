#BiocManager::install("GenomicRanges")
#BiocManager::install("genomation")
#install.packages("refGenome")
#install.packages(file.choose(), repos=NULL)
library(GenomicRanges)
library(rtracklayer)
library(devtools)
#library(refGenome)
#library(genomation)
library(dplyr)
#library(Gviz)
library(Biostrings)
library(ORFik)
library(optparse)
library(plyranges)

option_list = list(
  make_option(c("-a", "--assembly"), type="character", default=NULL, 
              help="Required -assembly name", metavar="assembly_name"),
  make_option(c("-c", "--chromosome"), type="character", default=NULL, 
              help="Required -chromosome number", metavar="chromosome")		  )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

File_Name = opt$assembly
Chomosome_number=opt$chromosome
Base_Directory=getwd()

###Starting read at Annotation_Folder
Folder_Name=paste0(Base_Directory,"/",File_Name,"/Annotation_steps")
setwd(Folder_Name)
getwd()
Masked_Fasta_File <- paste0(File_Name,"_",Chomosome_number,".masked.fasta")
Gff3_file <- paste0(File_Name,".gff3")

#STEP 1 mask the genome with bedtools maskfasta -fi FASTA -bed GFF -fo output
my_columns <- c("seqid","start","end","strand","type")
my_filter <- list(type="gene")

###Import GFF3 Genes HAN_Assembly and convert to Granges Object
gff.HAN_Assembly_raw_genes <- readGFF(Gff3_file,columns=my_columns, filter = my_filter)

Table.HAN_Assembly <- as.data.frame(gff.HAN_Assembly_raw_genes)
Table.HAN_Assembly$Chr <- Table.HAN_Assembly$seqid

HAN_Assembly_gr <- with(Table.HAN_Assembly, GRanges(Chr, IRanges(start, end),strand, id = Name))
HAN_Assembly_gr <- sortSeqlevels(HAN_Assembly_gr)
HAN_Assembly_gr <- sort(HAN_Assembly_gr, ignore.strand = TRUE)

###Starting read at EDTA folder
Folder_Name=paste0(Base_Directory,"/",File_Name,"/EDTA_Files")
setwd(Folder_Name)

###Import Fasta file and transform into DNAstring Dataset,names correlated with GFF3
Assembly_HAN_Assembly <- readDNAStringSet(Masked_Fasta_File)

Assembly_HAN_Assembly_names <- as.data.frame(Assembly_HAN_Assembly@ranges@NAMES)
colnames(Assembly_HAN_Assembly_names) <- "seqid"
Assembly_HAN_Assembly_names$Chr <- sub(" .*", "", Assembly_HAN_Assembly_names$seqid)

Assembly_HAN_Assembly@ranges@NAMES <- Assembly_HAN_Assembly_names$Chr

###Filter in order to split gff3 file

counter = as.data.frame(names(Assembly_HAN_Assembly))
colnames(counter) = "Chr"
Chr_HAN_Assembly_gr.table <- as.data.frame(HAN_Assembly_gr)
Chr_HAN_Assembly_gr.table$Chr <- Chr_HAN_Assembly_gr.table$seqnames

Chr_HAN_Assembly_gr.table_merged <- merge(counter,Chr_HAN_Assembly_gr.table, all.x=TRUE, by ="Chr")
Chr_HAN_Assembly_gr.table_merged <- filter(Chr_HAN_Assembly_gr.table_merged,is.na(id)==FALSE)
Chr_HAN_Assembly_gr<- with(Chr_HAN_Assembly_gr.table_merged, GRanges(Chr, IRanges(start, end),strand, id = id))

All_genes <- getSeq(Assembly_HAN_Assembly,Chr_HAN_Assembly_gr)

Masked_genes <- as.data.frame(alphabetFrequency(All_genes))
Masked_genes$length <- seqlengths(All_genes)
Masked_genes$TE_Overlap_Percentage <- (Masked_genes$N/Masked_genes$length)*100
Masked_genes$Name <- Chr_HAN_Assembly_gr@elementMetadata@listData$id

#Filter Genes less than 100bp and overlap below 75%
Masked_genes_filtered <- filter(Masked_genes, TE_Overlap_Percentage <=74.99)
Masked_genes_filtered <- filter(Masked_genes_filtered, length>=150)
Masked_genes_filtered_list <- Masked_genes_filtered[,c("Name","TE_Overlap_Percentage")]

HAN_Assembly_Curated_file <- merge(gff.HAN_Assembly_raw_genes,Masked_genes_filtered_list, by = "Name", all.y = TRUE)

###Starting read at Annotation_Folder
Folder_Name=paste0(Base_Directory,"/",File_Name,"/Annotation_steps")
setwd(Folder_Name)
getwd()
Output_name = paste0(File_Name,"_step1_output.gff3")
print(Output_name)

export(HAN_Assembly_Curated_file,Output_name,format="gff3")
