#/home/jmlazaro/scratch/Sunflower_annotation_Snakemake/HanHA89/Summary_data/HanHA89_busco/run_eudicots_odb10/missing_busco_list.tsv
#install.packages("tidyverse")
library(tidyverse)
library(jsonlite)
library(data.table)


busco_ids <- c("101979at71240","107327at71240")
names(busco_ids) <- busco_ids # so map_df gets .id right

# a tibble of interpro profile associated with each busco_id
busco_ipr <- map_df(busco_ids, .id="busco_id",  function(busco_id){
  write(busco_id, stderr()) # just so we can monitor progress
  # map the BUSCO ID to OrthoDB group ID
  query <- read_json(paste0("https://www.orthodb.org/search?query=", busco_id))
  odb_id <- query$data[1]
  # get all info on the orthogroup
  odb_info <- read_json(paste0("https://www.orthodb.org/group?id=", odb_id),
                        simplifyVector = TRUE)
  # return the FASTA of the missing BUSCO
  odb_fasta_info <- as.vector(unlist(fread(paste0("https://www.orthodb.org/fasta?id=", odb_id))))

  #Retunr the value into the function
  odb_info$data$interpro_domains
})

busco_ipr %>% select(busco_id, description)

