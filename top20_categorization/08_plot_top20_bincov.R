# # dlu 28.2.2023
# 
# this is basically the same thing as in bincov_plots BUT we changed two things:
#     1. take the id=95 mapping
#     2. viral regions use checkv in order to make the viral regions smaller (and more precise)
#     
# The script that created the old plots is backed up in ~/bioinf/


# author: dlueckin
# date: Tue Feb 28 12:12:50 2023

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


# import top20 ------------------------------------------------------------

top20_df <- fread("top20/top20_df_automatic_label.tsv")


# plot_function -----------------------------------------------------------

plotCov <- function(MAG, station, label){
    bincov_file <-  paste0("top20/detailed_MAG_mappings_95/bincov/", station, "_vs_", MAG, ".bincov")
    df <- fread(bincov_file, skip = 2)
    
    # viral plots are different
    if(label == "viral" | label == "contradicting"){
        checkv_file <-  paste0("top20/checkv/", MAG, "/contamination.tsv")
        
        checkv_df <- fread(checkv_file) %>% 
            filter(str_detect(string = region_types, pattern = "viral"))
        
        # if this is empty, create fake boundaries_df
        if(nrow(checkv_df) < 1){
            boundaries_df <- data.table(start = as.numeric(),
                                        end = as.numeric())
            label <- paste0(label, " - but no virus left after filtering")
            
        }else{
            # else we need do to the magic and get new region coordinates
            # create the boundaries df
            boundaries_df <- data.table()
            
            # go over the checkv_df, extract each region
            for(j in 1:nrow(checkv_df)){
                tmp_df <- data.table(big_region = checkv_df$contig_id[j],
                                     regions = unlist(str_split(checkv_df$region_types[j], pattern = ",")),
                                     boundaries = unlist(str_split(checkv_df$region_coords_bp[j], pattern = ",")),
                                     rel_start = 0,
                                     rel_end = 0)
                
                tmp_df$rel_start <- as.numeric(str_remove(tmp_df$boundaries, "-.*$"))
                tmp_df$rel_end <- as.numeric(str_remove(tmp_df$boundaries, "^\\d*-"))
                
                boundaries_df <- rbind(boundaries_df, tmp_df)
            }
            rm(tmp_df)
            
            # keep only viral regions
            boundaries_df <- boundaries_df %>% 
                filter(regions == "viral")
            
            # load old boundaries file, connect each big_region with original boundaries
            old_file <- paste0("top20/vs2/", MAG, "/final-viral-boundary.tsv")
            old_df <- fread(old_file) %>% 
                select(trim_bp_start, seqname_new)
            
            boundaries_df$big_region_start <- old_df$trim_bp_start[match(boundaries_df$big_region, old_df$seqname_new)]
            boundaries_df$start <- boundaries_df$big_region_start + boundaries_df$rel_start
            boundaries_df$end <- boundaries_df$big_region_start + boundaries_df$rel_end
            
            # NOW filter the regions by score, for that load score and connect
            score_file <- paste0("top20/vs2/", MAG, "/final-viral-score.tsv")
            score_df <- fread(score_file) %>% 
                select(seqname, max_score) %>% 
                filter(max_score >= 0.9)
            boundaries_df$score <- score_df$max_score[match(boundaries_df$big_region, score_df$seqname)]
            
            boundaries_df <- boundaries_df %>% 
                filter(score >= 0.9)
        }
    }else{
        # need fake boundaries_df
        boundaries_df <- data.table(start = as.numeric(),
                                    end = as.numeric())
    }
    
    # plot
    p <- ggplot(df, aes(x = RunningPos, y = Cov)) +
        geom_line() + 
        geom_segment(data = boundaries_df,
                     aes(x = start, xend = end,
                         y = 0, yend = 0),
                     color = "#FED766", size = 2) +
        ggtitle(paste0(station, "_vs_", MAG), subtitle = paste0("Mean Cov: ", round(mean(df$Cov), 2))) +
        theme_classic() +
        xlab("Position (bp)") +
        ylab("Coverage") +
        theme(legend.background = element_rect(fill = "transparent"),
              legend.box.background = element_rect(fill = "transparent"),
              panel.background = element_rect(fill = "transparent"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.background = element_rect(fill = "transparent", color = NA))
    # and return that plot
    return(p)
    
}


# plot and save -----------------------------------------------------------

for(i in 1:nrow(top20_df)){
    p <- plotCov(MAG = top20_df$name[i],
                 station = top20_df$station[i],
                 label = top20_df$automatic_label[i])
    
    ggsave(plot = p, filename = paste0("top20/bincov_plots_95/",
                                       top20_df$automatic_label[i], "/",
                                       top20_df$station[i], "_vs_", top20_df$name[i], ".png"),
           height = 5,
           width = 10,
           units = "cm")
}











# unimportant, just me trying to figrue something out ---------------------

df <- fread("top20/detailed_MAG_mappings_95/basecov/122_DCM_vs_GCA_002692925.1_ASM269292v1_genomic_concatenated.fna.basecov") %>% 
    filter(Pos > 5129567 & Pos < 5129895)
df2 <- fread("top20/detailed_MAG_mappings/basecov/122_DCM_vs_GCA_002692925.1_ASM269292v1_genomic_concatenated.fna.basecov")%>% 
    filter(Pos > 5129567 & Pos < 5129895)
df3 <- fread("test_16S_removal/GCA_002692925_no_16S_reads.basecov")%>% 
    filter(Pos > 5129567 & Pos < 5129895)
df4 <- fread("test_16S_removal/GCA_002692925_no_16S_reads_99.basecov")%>% 
    filter(Pos > 5129567 & Pos < 5129895)

station = "122_DCM"
MAG = "GCA_002692925"


p <- ggplot(df, aes(x = Pos, y = Coverage)) +
    geom_line() + 
    ggtitle(paste0(station, "_vs_", MAG), subtitle = "0.95") +
    theme_classic() +
    xlab("Position (bp)") +
    ylab("Coverage") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA))

p2 <- ggplot(df2, aes(x = Pos, y = Coverage)) +
    geom_line() + 
    ggtitle(paste0(station, "_vs_", MAG), subtitle = "0.90") +
    theme_classic() +
    xlab("Position (bp)") +
    ylab("Coverage") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA))

p3 <- ggplot(df3, aes(x = Pos, y = Coverage)) +
    geom_line() + 
    ggtitle(paste0(station, "_vs_", MAG), subtitle = "no16S - 0.95") +
    theme_classic() +
    xlab("Position (bp)") +
    ylab("Coverage") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA))


ggplot(df4, aes(x = Pos, y = Coverage)) +
    geom_line() + 
    ggtitle(paste0(station, "_vs_", MAG), subtitle = "no16S - 0.99") +
    theme_classic() +
    xlab("Position (bp)") +
    ylab("Coverage") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA))


a <- p/p2/p3/p4


