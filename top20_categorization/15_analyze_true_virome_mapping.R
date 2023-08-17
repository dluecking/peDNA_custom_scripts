# author: dlueckin
# date: Wed Mar  1 15:09:29 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))



# calc function -----------------------------------------------------------

calcPercMapped <- function(station, id){
    file1 <- paste0("true_virome_mappings/", station, "_vs_combined_MAGs_outm_", id, ".reads1.fq")
    file2 <- paste0("true_virome_mappings/", station, "_vs_combined_MAGs_outm_", id, ".reads2.fq")
    
    if(file.exists(file1)){
        r1 <- system(paste0("wc -l ", file1), intern = TRUE)
        r1 <- as.numeric(unlist(str_split(r1, "\\s"))[1]) 
        
        r2 <- system(paste0("wc -l ", file2), intern = TRUE)
        r2 <- as.numeric(unlist(str_split(r2, "\\s"))[1]) 
        
        perc_mapped <- (r1+r2) / 4 / 2000000 * 100

    }else{
        perc_mapped <- NA
    }
    return(perc_mapped)
}

# calc for each station ---------------------------------------------------

df <- fread("helper_files/list_of_stations.txt", header = FALSE)
names(df) <- c("station")

df$perc_mapped_76 <- unlist(sapply(df$station, calcPercMapped, id = 76))
df$perc_mapped_95 <- unlist(sapply(df$station, calcPercMapped, id = 95))


# fwrite(df %>% select(station, perc_mapped3), file = "perc_mapped_virome_reads_vs_MAG_per_station.tsv")

