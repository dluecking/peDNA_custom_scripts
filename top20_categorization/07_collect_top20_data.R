# author: dlueckin
# Fri Sep 16 11:45:48 2022 

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)


# working directory -------------------------------------------------------

setwd("/run/user/1000/gvfs/sftp:host=linux-desktop-1.mpi-bremen.de,user=dlueckin/home/dlueckin/projects/misc/three_stations")


# import top20 df ---------------------------------------------------------

top20_df <- fread("top20/top20_df_mapping.tsv")


# add viral prediction ----------------------------------------------------

vs2 <- fread("top20/combined_top20_final-viral-score.tsv")
names(vs2) <- c("seqname", "dsDNAphage", "ssDNA", "max_score", "max_score_group", "length", "hallmark", "viral", "cellular")
vs2 <- vs2 %>% 
    filter(max_score >= 0.9)

vs2$MAG <- str_remove(vs2$seqname, "\\|.*")

# fill in top20
top20_df$viral <- FALSE
for(i in 1:nrow(top20_df)){
    top20_df$viral[i] <- any(str_detect(string = vs2$MAG, pattern = top20_df$name[i]))
}
rm(vs2, i)


# virus activity old approach ---------------------------------------------

top20_df$virus_active_old <- FALSE

virus_activity_df <- data.table()

for(file in list.files("top20/virus_activity/", pattern = "rpkm")){
    
    tmp_rpkm <- fread(paste0("top20/virus_activity/", file), skip = 4)
    
    tmp_rpkm <- tmp_rpkm %>% 
        filter(RPKM >= 1)
    
    tmp_rpkm$sample <- str_remove(file, "\\_virus\\_activity\\_rpkm.txt")
    virus_activity_df <- rbind(virus_activity_df, tmp_rpkm)
}

# fill in
for(i in 1:nrow(virus_activity_df)){
    
    current_MAG <- str_remove(virus_activity_df$`#Name`[i], pattern = "\\|\\|.*")
    current_sample <-  virus_activity_df$sample[i]
    
    top20_df[top20_df$name == current_MAG & top20_df$station == current_sample]$virus_active_old <- TRUE
}
rm(tmp_rpkm, virus_activity_df, i, current_MAG, current_sample, file)



# virus activity new ------------------------------------------------------

top20_df$virus_active_new <- FALSE

isActive <- function(MAG, station){
    # read basecov in
    basecov_file <- paste0("top20/detailed_MAG_mappings/basecov/", station, "_vs_", MAG, ".basecov")
    basecov_df <- fread(basecov_file, skip = 2)
    names(basecov_df) <- c("mag", "pos", "cov")
    
    # read viral boundaries and score
    viral_boundary_df <- fread(paste0("top20/vs2/", MAG, "/final-viral-boundary.tsv"))
    viral_score_df <- fread(paste0("top20/vs2/", MAG, "/final-viral-score.tsv"))  %>% 
        filter(max_score >= 0.9) 
    
    # setup constants
    ONE_OR_MORE_REGIONS_ACTIVE <- FALSE
    BASECOV_COV_AVG <- mean(basecov_df$cov)
    BASECOV_COV_SD <- sd(basecov_df$cov)
    
    # check each region, if cov is 2xSD above mean
    for(region in viral_score_df$seqname){
        region_start <- viral_boundary_df$prox_bp_start[viral_boundary_df$seqname_new == region]
        region_end <- viral_boundary_df$prox_bp_end[viral_boundary_df$seqname_new == region]
        
        region_df <- basecov_df %>% 
            filter(pos >= region_start & pos <= region_end)
        
        REGION_AVG_COV <- mean(region_df$cov)
        
        if(REGION_AVG_COV >= 2 * BASECOV_COV_SD + BASECOV_COV_AVG){
            ONE_OR_MORE_REGIONS_ACTIVE <- TRUE
        }
    }
    return(ONE_OR_MORE_REGIONS_ACTIVE)
}


# apply to each mag with a virus
# this takes time!
# for(i in 1:nrow(top20_df)){
#     if(top20_df$viral[i] == TRUE){
#         top20_df$virus_active_new[i] <- isActive(MAG = top20_df$name[i], station = top20_df$station[i])
#         
#     }
# }




# add GTA info ------------------------------------------------------------

gta <- fread("top20/gtas/blastp_result.out")
names(gta) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

gta$MAG <- str_remove(gta$sseqid, "\\_\\d*$")
gta <- gta %>% 
    filter(pident > 50 & evalue < 10^-5)

top20_df$no_of_GTA_hits <- 0

for(i in 1:nrow(top20_df)){
    top20_df$no_of_GTA_hits[i] <- gta %>% 
        filter(MAG == top20_df$name[i]) %>%
        distinct(qseqid) %>% 
        nrow() %>% 
        unlist()
}

top20_df$no_of_GTA_hits <- as.numeric(top20_df$no_of_GTA_hits)
top20_df$gta <- FALSE
top20_df$gta[top20_df$no_of_GTA_hits >= 10] <- TRUE

rm(gta, i)


# add meta data to each MAG -----------------------------------------------

metadata <- fread("PRJNA391943_AssemblyDetails.txt")
names(metadata)  <- c("Assembly", "Level", "WGS", "BioSample", "Isolate", "Taxonomy", "NOTHING")

metadata <- metadata %>% 
    select(-NOTHING)


top20_df$Assembly <- str_remove(top20_df$name, "\\_ASM.*$")
top20_df <- left_join(top20_df, metadata, by = "Assembly")


# add rpkm and coverage info ----------------------------------------------

rpkm_df <- data.table()

for(i in 1:nrow(top20_df)){
    rpkm_file <- paste0("top20/detailed_MAG_mappings/rpkm/", top20_df$station[i], "_vs_", top20_df$name[i], ".rpkm")
    tmp_df <- fread(rpkm_file, skip = 4)
    tmp_df$station <- top20_df$station[i]
    
    rpkm_df <- rbind(rpkm_df, tmp_df)
}

top20_df$tmp_id <- paste0(top20_df$name, "-", top20_df$station)
rpkm_df$tmp_id <- paste0(rpkm_df$`#Name`, "-", rpkm_df$station)

rpkm_df <- rpkm_df %>% 
    select(tmp_id, Length, Coverage, Reads)

top20_df <- left_join(top20_df, rpkm_df, by = c("tmp_id"))



# final selection ---------------------------------------------------------
top20_df$automatic_label <- "ev_producer"

# if virus and some coverage
top20_df[viral == TRUE & gta == FALSE & virus_active_old == TRUE, "automatic_label" := "viral"]

# if virus and very high coverage
top20_df[viral == TRUE & gta == FALSE & virus_active_new == TRUE, "automatic_label" := "viral"]

# if gta > 10
top20_df[gta == TRUE, "automatic_label" := "gta"]

# if gta somehwere 5-9 and viral 
top20_df[viral == TRUE & no_of_GTA_hits > 5 & no_of_GTA_hits < 10,  "automatic_label" := "contradicting"]




# write -------------------------------------------------------------------

top20_df_to_write <- top20_df %>%
    select(name, station, automatic_label, gta, no_of_GTA_hits, viral, virus_active_new, virus_active_old, 
           Taxonomy, Length, Coverage, Reads, assignedReads, assignedBases) %>%
    arrange(automatic_label)

names(top20_df_to_write) <- c("name", "station", "automatic_label", "gta", "no_of_GTA_hits", "viral", "virus_active_new", "virus_active_old",
                              "taxonomy", "length", "coverage", "high_id_reads", "low_id_reads", "low_id_bases")
 
fwrite(top20_df_to_write, "top20/top20_df_automatic_label.tsv")




# # writing out lists and df ------------------------------------------------
# 
# 
# # fwrite(top20_df %>% filter(final_label == "gta"), "final_top20_df.tsv", append = TRUE)
# 
# viruses <- top20_df %>% group_by(sample) %>% filter(final_label == "viral") %>% arrange(by = desc(Coverage))
# # fwrite(viruses, "final_top20_df.tsv", append = TRUE)
# 
# evs <- top20_df %>% filter(final_label == "ev_producer") %>% group_by(sample) %>% slice_max(order_by = Coverage, n = 3)
# # fwrite(evs, "final_top20_df.tsv", append = TRUE)
# 
# 
# 
# 
# 
# ################################################################################
# # THIS PART WAS NEEDED, but is not we do not run it again ######################
# ################################################################################
# 
# 
# # # is the viral region active? ---------------------------------------------
# # h <- top20_df
# # h$`#Name` <- str_replace(h$`#Name`, "\\.fna", "\\_100kb\\_filtered\\.fna")
# # 
# # # I have data for 122_MES, 158_SRF and 64_DCM
# # to_write <- h %>% 
# #     filter(sample == "122_MES") %>% 
# #     filter(viral == TRUE) %>% 
# #     select(`#Name`) %>% 
# #     unlist()
# # write(to_write, file = "virus_activity/map_the_viral_seqs_of_these_MAGS_to_122_MES.txt")
# # 
# # 
# # to_write <- h %>% 
# #     filter(sample == "158_SRF") %>% 
# #     filter(viral == TRUE) %>% 
# #     select(`#Name`) %>% 
# #     unlist()
# # write(to_write, file = "virus_activity/map_the_viral_seqs_of_these_MAGS_to_158_SRF.txt")
# # 
# # 
# # to_write <- h %>% 
# #     filter(sample == "64_DCM") %>% 
# #     filter(viral == TRUE) %>% 
# #     select(`#Name`) %>% 
# #     unlist()
# # write(to_write, file = "virus_activity/map_the_viral_seqs_of_these_MAGS_to_64_DCM.txt")
# 






