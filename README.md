## Extracellular vesicles are the main contributor to the non-viral protected extracellular sequence space
## Supplementary Custom Scripts
Scripts, commands and programs used to the manuscript "Extracellular vesicles are the main contributor to the non-viral protected extracellular sequence space".
There are two cases where custom scripts were used and referenced as such in the manuscript.
1. Summarizing the results of a viral identification pipeline according to standards established elsewhere (https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3).
2. Categorization of MAGs depending on their main mechanism of horizontal gene transfer, a result of this study. The critical categorization script is `top20_categorization/07_collect_top20_data.R`, in which the results of previous steps are collected and MAGs are labelled according to the following logic:
Excerpt of `07_collect_top20_data.R`:
```R
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
```
Keep in mind, that this resulted in a 'automatic label' which then was manually curated by scrutinizing coverage plots, verifying the VS2 viral prediction score. 'Contradicting' cases were resolved by verifying the activity (read: increased coverage of predicted viral region) of a given phage. 


