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

stats_files <- snakemake@input[["stats_files"]]
plot_file <- snakemake@output[["plot"]]

# dev
# stats_files <- list.files("output/050_stats",
#            recursive = TRUE,
#            pattern = "stats.tsv",
#            full.names = TRUE)
# plot_file <- "test/stats.pdf"

########
# MAIN #
########

# read stats
names(stats_files) <- gsub("^output/050_stats/(.+)/stats.tsv", "\\1", stats_files)
stats_list <- lapply(stats_files, fread)
stats <- rbindlist(stats_list, idcol = "filename")

# mung
stats[, mp := dirname(filename)]
stats[, c("read_set", "k", "diplo") := tstrsplit(basename(filename), "_")]

# fill in missing values during melt. key order has to match order in CJ()
setkey(stats, mp, read_set, k, diplo)
pd <- melt(stats[CJ(unique(mp), unique(read_set), unique(k), unique(diplo))],
           id.vars = c("mp", "read_set", "k", "diplo"),
           measure.vars = c("n_scaffolds", "scaf_bp", "scaf_N50", "scaf_L50"),
           fill = TRUE)

# set up labels
read_set_order <- c( "raw" = "Raw",
                     "norm" = "Normalised")
variable_order <- c("n_scaffolds" = "Scaffolds (K)",
                    "scaf_bp" = "Assembled size (Mb)",
                    "scaf_N50" = "N50 (K)",
                    "scaf_L50" = "L50 (Kb)")

pd[, value := as.double(value)]
pd[variable == "n_scaffolds", value := value / 1000]
pd[variable == "scaf_bp", value := value / 1e6]
pd[variable == "scaf_N50", value := value / 1000]
pd[variable == "scaf_L50", value := value / 1000]


pd[, k := as.numeric(gsub("[^[:digit:]]+", "", k))]
pd[, diplo := as.numeric(gsub("[^[:digit:]]+", "", diplo))]

pd[, read_set := factor(plyr::revalue(read_set, read_set_order),
                        levels = read_set_order)]

pd[, variable := factor(plyr::revalue(variable, variable_order),
                        levels = variable_order)]

# set up colours
fill_colours = viridis::viridis_pal()(3)

# line for max per value
best_dt <- pd[, .(value = ifelse(
    variable %in% c("Scaffolds (K)", "N50 (K)"),
    min(value, na.rm = TRUE),
    max(value, na.rm = TRUE))),
    by = variable]

# draw the plot
gp <- ggplot(pd, aes(x = as.factor(k),
                     y = value,
                     fill = as.factor(diplo))) +
    coord_flip() +
    theme(strip.placement = "outside", strip.background.x = element_blank()) +
    facet_grid(mp + read_set ~ variable, scales = "free", switch = "x") +
    xlab(expression(italic("k"))) + ylab(NULL) +
    scale_fill_manual(values = fill_colours[1:2],
                      guide = guide_legend(title = "Ploidy")) +
    geom_col(position = position_dodge(width = 0.8),
             width = 0.5) +
    geom_hline(mapping = aes(yintercept = value),
               data = best_dt,
               colour = fill_colours[3])

# write output
ggsave(filename = plot_file,
       plot = gp,
       device = cairo_pdf,
       width = 10,
       height = 7.5,
       units = "in")

# write session info
sessionInfo()
