# author: dlueckin
# Wed Jul 13 12:41:20 2022 

  
# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)


# working directory -------------------------------------------------------
setwd("/run/user/1000/gvfs/sftp:host=linux-desktop-1.mpi-bremen.de,user=dlueckin/home/dlueckin/projects/misc/three_stations/")


# import data -------------------------------------------------------------

stations <- fread("helper_files/list_of_stations.txt", header = F)

data <- data.table()

for(station in stations$V1){
    tmp_df <- rbindlist(lapply(list.files("mappings",
                                          pattern = station,
                                          full.names = TRUE),
                               fread))
    tmp_df$station <- station
    data <- rbind(data, tmp_df)
}
rm(tmp_df)

names(data) <- c("name", "%unambiguousReads", "unambiguousMB", "%ambiguousReads",
                 "ambiguousMB", "unambiguousReads", "ambiguousReads",
                 "assignedReads", "assignedBases", "station")





# select_top20_each -------------------------------------------------------

a <- fread("../../mvome_pipeline/manuscript_plots/input-data - top20_combined.csv") %>% 
    filter(Sample != "heligoland")

top20_df <- data %>% 
    group_by(station) %>% 
    slice_max(order_by = assignedReads, n = 20)

fwrite(top20_df, "top20/top20_df_mapping.tsv")
write(top20_df$name, file = "helper_files/top20_contigs.txt")
write(unique(top20_df$name), file = "helper_files/top20_contigs_unique.txt")


# for(i in 1:nrow(a)){
#     n <- a$shortname[i]
#     if(any(str_detect(top20_df$name, n))){
#         print("found!")
#     }else{
#         print("not found:")
#         print(n)
#     }
# }





