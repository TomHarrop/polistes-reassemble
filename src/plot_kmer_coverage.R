#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")


library(data.table)
library(bit64)
library(ggplot2)
library(scales)

#############
# FUNCTIONS #
#############

GenerateHistfileName <- function(histfile) {
    my_file <- sub(".txt", "", basename(histfile))
    my_hist <- ifelse(grepl("_out$", my_file), "Normalised", "Raw")
    my_lib <- gsub("-.*", "", my_file)
    return(paste(my_lib, my_hist, sep = "-"))}


###########
# GLOBALS #
###########

hist_before_file <- snakemake@input[["hist_before"]]
hist_after_file <- snakemake@input[["hist_after"]]
peak_file <- snakemake@input[["peaks"]]
plot_file <- snakemake@output[["plot"]]

# find files
hist_files <- list.files("output/020_norm",
                         pattern = "hist.*txt",
                         full.names = TRUE)
peak_file <- "output/020_norm/500_bp_insert_standard-peaks.txt"

# dev
# hist_before_file <- "output/030_norm/hist.txt"
# hist_after_file <- "output/030_norm/hist_out.txt"
# peak_file <- "output/030_norm/peaks.txt"

########
# MAIN #
########

# read data
peaks <- fread(paste("grep '^[^#]'", peak_file))
names(hist_files) <- GenerateHistfileName(hist_files)
hist_data_list <- lapply(hist_files, fread)
combined_data <- rbindlist(hist_data_list, idcol = "histfile_name")
combined_data[, c("lib", "type") := tstrsplit(histfile_name, "-", fixed = TRUE)]

# arrange plot
all_libs <- rev(combined_data[, unique(lib)])
lib_order <- structure(gsub("_", " ", all_libs), names = all_libs)
combined_data[, lib := factor(plyr::revalue(lib, lib_order),
                              levels = lib_order)]
combined_data[, type := factor(type, levels = c("Raw", "Normalised"))]

# hlines
mincov <- peaks[1, V1]
p1 <- peaks[1, V2]
maxcov <- peaks[1, V3]

# plot title
gt <- paste0(
    p1, "× 31-mer coverage. ",
    "Main peak: ", mincov, "×–", maxcov, "×"
)

# plot
kmer_plot <- ggplot(combined_data, aes(x = `#Depth`,
                                       y = Unique_Kmers,
                                       colour = lib,
                                       linetype = type)) +
    theme_minimal() +
    theme(legend.position = c(5/6, 2/4)) +
    geom_vline(xintercept = c(mincov, p1, maxcov),
               colour = "grey50") +
    geom_path(alpha = 0.75) +
    scale_colour_viridis_d(guide = guide_legend(title = NULL)) +
    scale_linetype_discrete(guide = guide_legend(title = NULL)) +
    scale_y_continuous(
        trans = "log10",
        labels = trans_format("log10", math_format(10^.x)),
        breaks = trans_breaks("log10", function(x) 10^x)) +
    scale_x_continuous(trans = log_trans(base = 4),
                       breaks = trans_breaks(function(x) log(x, 4),
                                             function(x) 4^x)) +
    xlab("31-mer depth") + ylab("Number of unique 31-mers") + ggtitle(gt)

# write output
ggsave(filename = plot_file,
       plot = kmer_plot,
       device = cairo_pdf,
       width = 10,
       height = 7.5,
       units = "in")

# write session info
sessionInfo()
