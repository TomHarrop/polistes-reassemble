#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(ggplot2)

###########
# GLOBALS #
###########


busco_files <- snakemake@input[["busco_files"]]
plot_file <- snakemake@output[["plot"]]

# dev
# plot_file <- "output/050_assembly-stats/assembly_stats.pdf"
# busco_files <- list.files("output/060_busco",
#                           recursive = TRUE,
#                           pattern = "full_table_busco.tsv",
#                           full.names = TRUE)

########
# MAIN #
########

# read data
names(busco_files) <- gsub("^output/060_busco/(.+)/run_busco/full_table_busco.tsv",
                           "\\1",
                           busco_files)

busco_list <- lapply(busco_files, fread, fill = TRUE, skip = 4)
busco_data <- rbindlist(busco_list, idcol = "filename")

# mung
busco_data[, c("mp", "read_set", "k", "diplo") := tstrsplit(filename, "/")]

# summarise
busco_counts <- busco_data[, .(busco_count = length(unique(`# Busco id`))),
                           by = .(mp, read_set, k, diplo, Status)]
busco_counts[, total_buscos := sum(busco_count), by = .(mp, read_set, k, diplo)]
busco_counts[, busco_percent := busco_count * 100 / total_buscos]

# fill blanks
setkey(busco_counts, mp, read_set, k, diplo, Status)
pd <- busco_counts[CJ(unique(mp),
                      unique(read_set),
                      unique(k),
                      unique(diplo),
                      unique(Status))]

# order the plot
read_set_order <- c( "raw" = "Raw",
                     "norm" = "Normalised")
status_order <- c("Complete", "Duplicated", "Fragmented", "Missing")

pd[, k := as.numeric(gsub("[^[:digit:]]+", "", k))]
pd[, diplo := as.numeric(gsub("[^[:digit:]]+", "", diplo))]

pd[, read_set := factor(plyr::revalue(read_set, read_set_order),
                        levels = read_set_order)]

pd[, Status := factor(Status, levels = rev(status_order))]

# set up colours
fill_colours = viridis::viridis_pal()(5)

# winning assembly completeness
max_pct <- pd[Status == "Complete", max(busco_percent, na.rm = TRUE)]

# plot
gp <- ggplot(pd, aes(x = as.factor(k),
                     y = busco_percent,
                     fill = Status)) +
    coord_flip() +
    ylab("%") + xlab("k") +
    scale_fill_manual(values = rev(fill_colours[2:5])) +
    facet_grid(mp + read_set ~ diplo) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = max_pct, colour = fill_colours[1])

# write output
ggsave(filename = plot_file,
       plot = gp,
       device = cairo_pdf,
       width = 10,
       height = 7.5,
       units = "in")

# log
sessionInfo()
