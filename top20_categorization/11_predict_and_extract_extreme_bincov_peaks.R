# author: dlueckin
# date: Tue Jan 17 11:53:45 2023

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


# predict function --------------------------------------------------------

for(f in list.files("top20/detailed_MAG_mappings_95/bincov/", full.names = TRUE)){
    df <- fread(f, skip = 2)
    peak <- which(which(df$Cov >= (mean(df$Cov) + 20*sd(df$Cov))))
    
    
    
    
    
}
