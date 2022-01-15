#bedtools maskfasta -fi ../Ref/PK_Scaffold_Split/PK_hap2.scaffold_10.fasta -bed PK_hap2.scaffold_10.fasta.mod.EDTA.intact.gff3 -fo PK_hap2.scaffo ld_10.masked.fasta              

library(GenomicFeatures)
library(Biostrings)
library(dplyr)
library(ORFik)
library(optparse)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome)

#usage: R --vanilla < Step2_Filtering.R
#--args -g HAN_Filtered_PK_hap2.scaffold_10.all.AED_0.8.gff3 -a PK_hap2.scaffold_10.fasta -m PK_hap2.scaffold_10.masked.fasta -o PK_hap2.scaffold_10

#https://cran.r-project.org/web/packages/optparse/readme/README.html

option_list = list(
  make_option(c("-g", "--gff_file"), type="character", default=NULL, 
              help="Required -GFF3 file", metavar="Fasta_File"),
  make_option(c("-a", "--assembly"), type="character", default=NULL, 
              help="Required unmasked Assembly Name file", metavar="Assembly_File"),
  make_option(c("-m", "--massembly"), type="character", default=NULL, 
              help="Required masked Assembly Name file", metavar="Assembly_File"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Output Name file", metavar="Assembly_File"),
  make_option(c("-s", "--start_gff3"), type="character", default=NULL, 
              help="Starting GFF3", metavar="Starting_GFF3_File")			  
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

File_Name = opt$gff_file
Assembly_Name = opt$assembly
Masked_Assembly_Name = opt$massembly
Output_Name = opt$output
Original_GFF3 = opt$start_gff3

my_columns <- c("seqid","start","end","strand","type")
my_filter <- list(type="gene")

HAN_scaffold_gff3_Original <- as.data.frame(readGFF(Original_GFF3))
		#saveRDS(HAN_scaffold_gff3_Original, file = "HAN_scaffold_gff3_Original.rds")

HAN_scaffold_gff3_genes<- as.data.frame(readGFF(File_Name,columns = my_columns, filter = my_filter))
		#saveRDS(HAN_scaffold_gff3_genes, file = "HAN_scaffold_gff3_genes.rds")	

###Import Fasta file and transform into DNAstring Dataset,names correlated with GFF3
Assembly_DNA <- readDNAStringSet(Assembly_Name)
		#saveRDS(Assembly_DNA, file = "Assembly_DNA.rds")
		
Masked_Assembly_DNA <- readDNAStringSet(Masked_Assembly_Name)
		#saveRDS(Masked_Assembly_DNA, file = "Masked_Assembly_DNA.rds")

HAN_Genes_gr <- with(HAN_scaffold_gff3_genes, GRanges(seqid, IRanges(start, end),strand, id = Name))



HAN_All_genes <- getSeq(Masked_Assembly_DNA,HAN_Genes_gr)
HAN_Masked_genes <- as.data.frame(alphabetFrequency(HAN_All_genes))
HAN_Masked_genes$length <- seqlengths(HAN_All_genes)
HAN_Masked_genes$TE_Overlap_Percentage <- (HAN_Masked_genes$N/HAN_Masked_genes$length)*100
HAN_Masked_genes$Name <- HAN_Genes_gr@elementMetadata@listData$id


#Filter Genes less than 100bp and overlap below 75%
HAN_Masked_genes_filtered <- filter(HAN_Masked_genes, TE_Overlap_Percentage <=74.99)
HAN_Masked_genes_filtered <- filter(HAN_Masked_genes_filtered, length>=150)
HAN_Masked_genes_filtered_list <- HAN_Masked_genes_filtered[,c("Name","TE_Overlap_Percentage")]

HAN_genes_Curated_file <- merge(HAN_scaffold_gff3_genes,HAN_Masked_genes_filtered_list, by = "Name", all.y = TRUE)
####END STEP1

HAN_Original.gr <- makeGRangesFromDataFrame(HAN_scaffold_gff3_Original, ignore.strand = FALSE,
                                          seqinfo = NULL,
                                          seqnames.field = "seqid",
                                          start.field = "start",
                                          end.field = "end",
                                          strand.field = "strand",
                                          keep.extra.columns = TRUE)

										  
step1_filtered_genes.gr <- makeGRangesFromDataFrame(HAN_genes_Curated_file , ignore.strand = FALSE,
                                                    seqinfo = NULL,
                                                    seqnames.field = "seqid",
                                                    start.field = "start",
                                                    end.field = "end",
                                                    strand.field = "strand",
                                                    keep.extra.columns = TRUE)

HAN_Filtered_GFF <- subsetByOverlaps(HAN_Original.gr, step1_filtered_genes.gr)

###CDS Filtering
gff.HAN_CDS <- filter(as.data.frame(HAN_Filtered_GFF),type=="CDS")
HAN_CDS_gr <- with(gff.HAN_CDS, GRanges(seqnames, IRanges(start, end),strand, id = ID))


#ORF Analysis START
##################
## make TxDb from GTF file 
txdb <- makeTxDbFromGRanges(HAN_Filtered_GFF)

## get intron information
all.introns <- intronicParts(txdb)
	#saveRDS(all.introns, file = "all.introns.rds")
all.exons <- exonicParts(txdb)
	#saveRDS(all.exons, file = "all.exons.rds")

ORF_filtered_genes.gr <-GRanges()
Proteins_list <- list()

for (i in 1:length(step1_filtered_genes.gr)) 
{
  Gene_Number <- i
  print(Gene_Number)
  Sequence_to_Analyze <- getSeq(Masked_Assembly_DNA, step1_filtered_genes.gr[Gene_Number])
  Sequence_to_Analyze <- DNAString(as.character(Sequence_to_Analyze))
  
  Sequence_Strand <- as.character(step1_filtered_genes.gr[Gene_Number]@strand)
  Reverse_CDS=FALSE
  
  #Reverse complement if in the negative Strand
  if (Sequence_Strand=="-"){
    Sequence_to_Analyze <- reverseComplement(Sequence_to_Analyze)
    Reverse_CDS=TRUE
  }
  
  #ORF bigger than 150bp (50aa)
  Possible_ORFS <- as.data.frame(findORFs(as.character(Sequence_to_Analyze)))
  
  Possible_ORFS <- filter(Possible_ORFS, width > 150)
  
  Gene_Name <- as.character(step1_filtered_genes.gr[Gene_Number]@seqnames)
  Gene_Pos <- step1_filtered_genes.gr[Gene_Number]@ranges
  #ORFS count
  No_ORFS <- as.numeric(nrow(Possible_ORFS))
  
  #Exons and Introns Count
  Lookup_Gene <- as.data.frame(step1_filtered_genes.gr[Gene_Number])
  
  Lookup_Gene <- with(Lookup_Gene, GRanges(seqnames, IRanges(start, end),strand, id))
  
  #table(!is.na(findOverlaps(all.exons, Lookup_Gene, select="arbitrary")))
  
  CDS_region <- as.data.frame(subsetByOverlaps(HAN_CDS_gr, Lookup_Gene))
  
  CDS_region <- with(CDS_region, GRanges(seqnames, IRanges(start, end),strand, id))
  ###Order the SEQUENCE
  CDS_region <- sortSeqlevels(CDS_region)
  CDS_region <- sort(CDS_region, decreasing = Reverse_CDS)
      ## Drop all unused seqlevels:
      seqlevels(CDS_region) <- seqlevelsInUse(CDS_region)
  
  if(length(CDS_region)>0)
  {
    CDS_Sequence <- getSeq(Assembly_DNA,CDS_region)
    CDS_Sequence <- do.call(xscat, CDS_Sequence)
    
    Protein_translated <- as.character(translate(CDS_Sequence, if.fuzzy.codon="solve"))
  } else
  {
    Protein_translated <- "NNNNNNNNNNNNNNNNNNNNNNN"
  }
  #Adding to protein list
  Proteins_list <- append(Proteins_list,c(paste0(">",step1_filtered_genes.gr[Gene_Number]$ID),Protein_translated))
  
  Exons_Count <- as.numeric(nrow(as.data.frame(subsetByOverlaps(all.exons, Lookup_Gene))))
  
  #table(!is.na(findOverlaps(all.introns, Lookup_Gene, select="arbitrary")))
  Introns_Count <- as.numeric(nrow(as.data.frame(subsetByOverlaps(all.introns, Lookup_Gene))))
  
  #Granges Object Creations
  Gene_Granges <- GRanges(seqnames = Gene_Name,ranges = Gene_Pos, strand = Sequence_Strand, 
                          ORFS_Number = No_ORFS, No_Exons = Exons_Count, No_Introns = Introns_Count )
  ORF_filtered_genes.gr <- append(ORF_filtered_genes.gr,Gene_Granges)
}

ORF_filtered_genes.DF <- as.data.frame(ORF_filtered_genes.gr)
Pseudogenes <- ORF_filtered_genes.DF

ORF_filtered_genes.DF <- filter(ORF_filtered_genes.DF,ORFS_Number >0)

ORF_filtered_genes.gr <- with(ORF_filtered_genes.DF, GRanges(seqnames, IRanges(start, end),
                                                             strand, ORFS_Number))

###Getting Step2 Sequence filtered GFF3
HAN_Step2_Filtered_GFF <- as.data.frame(subsetByOverlaps(HAN_Original.gr, ORF_filtered_genes.gr))
print(colnames(HAN_Step2_Filtered_GFF))
#HAN_Step2_Filtered_GFF <- HAN_Step2_Filtered_GFF[,c("seqnames","start","end","width","strand","source","type","score","phase","ID","Name","locus_tag")]

export(HAN_Step2_Filtered_GFF,paste0(Output_Name,".gff3"),format="gff3")

sink(paste0("predicted_proteins_HAN_ORF_Filtered_",Output_Name,".fasta"))
writeLines(unlist(lapply(Proteins_list, paste, collapse=" ")))
sink()

####Pseudogenes
Pseudogenes$Pseudo <- ifelse(Pseudogenes$ORFS_Number==0,"PseudoGene","Gene")
Pseudogenes$Pseudo <- ifelse(Pseudogenes$No_Introns==0,"Likely_PseudoGene",Pseudogenes$Pseudo)
Pseudogenes$Gene_Name <- step1_filtered_genes.gr$ID

export(HAN_Step2_Filtered_GFF,paste0("predicted_Pseudogenes_HAN_ORF_Filtered_",Output_Name,".gff3"),format="gff3")
