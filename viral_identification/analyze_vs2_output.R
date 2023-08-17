# title:
# author: domi
# date: Wed Jun  9 13:19:33 2021

# libraries ---------------------------------------------------------------
library(data.table)
suppressMessages(library(dplyr))

# import data -------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)

if(length(args) != 3){
    print("Usage:\nRscript analyze_vs2_output.R final-viral-score.tsv contamination.tsv <output_path>")
    quit()
}


vs_dt_path = args[1]
cv_dt_path = args[2]
output_path = args[3]

vs_dt <- fread(vs_dt_path)
cv_dt <- fread(cv_dt_path)

names(vs_dt)[1] <- "contig_id"
dt <- left_join(vs_dt, cv_dt, by = "contig_id")

rm(vs_dt, cv_dt)


# assign label ------------------------------------------------------------
# based on: https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-btv8nn9w?step=4

dt$label <- ""
dt$max_score[dt$max_score == "NaN"] <- 0
for(i in 1:nrow(dt)){
    if(dt$viral_genes[i] > 0 || dt$max_score[i] > 0.95 ||  dt$hallmark[i] > 2 || (dt$viral_genes[i] == 0 && dt$host_genes[i] == 0))
        dt$label[i] <- "true viral"

    if((dt$viral_genes[i] == 0 && dt$host_genes[i] > 1) || (dt$viral_genes[i] == 0 && dt$host_genes[i] == 1 && dt$length[i] < 10000))
        dt$label[i] <- "not viral"

    if((dt$viral_genes[i] == 0 && dt$host_genes[i] == 1 && dt$length[i] > 10000))
        dt$label[i] <- "manual curation"
}


# save dt -----------------------------------------------------------------
fwrite(dt, output_path, col.names = TRUE, sep = "\t")


# print report ------------------------------------------------------------

print("Quick Summary")
print(table(dt$label))
