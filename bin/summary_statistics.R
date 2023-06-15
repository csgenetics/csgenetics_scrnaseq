#! /usr/bin/env Rscript

# load libraries
library(Matrix, quietly=T)
library(rjson, quietly=T)
library(optparse)

# ------------------------------ Parse Arguements

arg_list <- list(
  make_option(
    c('--sample_id'), action = 'store', type = 'character',
    help = 'sample id'),
  make_option(
    c('--matrix'), action = 'store', type = 'character',
    help = 'matrix'),
  make_option(
    c('--barcodes'), action = 'store', type = 'character',
    help = 'barcode'),
  make_option(
    c('--features'), action = 'store', type = 'character',
    help = 'feature'),
  make_option(
    c('--multiqc_json'), action = 'store', type = 'character',
    help = 'Path to the multiqc json file. REQUIRED'),
  make_option(
    c('--antisense'), action = 'store', type = 'character',
    help = 'antisense argument'),
  make_option(
    c('--qualimap_report'), action = 'store', type = 'character',
    help = 'Path to the qualimap report txt  file. REQUIRED'),
  make_option(
    c('--cell_caller'), action = 'store', type = 'integer',
    help = 'The minimum number of nuclear genes detected for a barcode to be considered a cell. REQUIRED')
  )

opt <- parse_args(OptionParser(option_list=arg_list))

# ------------------------------ Load datasets

# load filter mtx
filtered_mtx <- readMM(opt$matrix) # load filtered mtx

# Read in filter `features.tsv`
features <- read.csv(opt$features, sep = '\t', header = F)
rownames(filtered_mtx) <- features[,2] # attach gene_ids

# Read in filter `barcodes.tsv`
cell_ids <- read.csv(opt$barcodes, sep = '\t', header = F)
colnames(filtered_mtx) <- cell_ids$V1

### load multiqc json output
data <- fromJSON(file=opt$multiqc_json)
total_number_reads   <- data$report_general_stats_data[[4]][[paste0(opt$sample_id,"_R1")]]$total_sequences

# ------------------------------ Table formatting Mapping and Seq statistics

valid_barcode <- data$report_general_stats_data[[3]][[paste0(opt$sample_id,"_R1")]]$before_filtering_total_reads

seq_sat <- paste0(round(data$report_general_stats_data[[3]][[paste0(opt$sample_id,"_R1")]]$pct_duplication,1),"%")

q30_read <- paste0(round(data$report_saved_raw_data$multiqc_fastp[[paste0(opt$sample_id,"_R1")]]$summary$before_filtering$q30_rate * 100,1),"%")

q30_barcode= paste0(round(data$report_saved_raw_data$multiqc_fastp[[paste0(opt$sample_id,"_R2")]]$summary$after_filtering$q30_rate * 100,1),"%")

### read_mapped_genome
mapped = data$report_general_stats_data[[2]][[opt$sample_id]]

read_mapped_genome =  paste0(mapped$uniquely_mapped+ mapped$multimapped," ","("
                            ,mapped$uniquely_mapped_percent+ mapped$multimapped_percent,"%",")")

read_uniq_mapped_genome = paste0(mapped$uniquely_mapped," ","(",mapped$uniquely_mapped_percent,"%",")")

quali <- readLines(opt$qualimap_report)

exonic <- sub(".*=\\s*(.*)", "\\1", grep("exonic =", quali, value = TRUE, fixed = TRUE))
intronic <- sub(".*=\\s*(.*)", "\\1", grep("intronic =", quali, value = TRUE, fixed = TRUE))
intergenic <- sub(".*=\\s*(.*)", "\\1", grep("intergenic =", quali, value = TRUE, fixed = TRUE))
totalAlignments <- sub(".*=\\s*(.*)", "\\1", grep("total alignments =", quali, value = TRUE, fixed = TRUE))

val_antisense <- read.csv(opt$antisense,sep = '\t', header = F)
per_antisense <- round((as.numeric(val_antisense$V1)/ as.numeric(gsub(",", "", totalAlignments))) * 100,2)
antisense <- paste0(val_antisense$V1," ","(",per_antisense,"%",")")

# ------------------------------ Write table Mapping and Sequencing statistics

seq_stats <- data.frame(stat = c('Total Number of Reads', 'Valid Barcodes / CSGX Cell Barcodes', 'Sequencing Saturation', # get sequencing/alignment stats 
                                  'Q30 Bases in Barcode','Q30 Bases in RNA Read'), 
                        value = prettyNum(c(total_number_reads, valid_barcode, seq_sat, q30_barcode,q30_read), big.mark = ','))
write.table(seq_stats, paste0(opt$sample_id,"_seq_stats.csv"),quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")

map_stats <- data.frame(stat = c('Reads Mapped to Genome', 'Reads Mapped Confidently to Genome', 'Reads Mapped Confidently to Intergenic Regions',  
                                 'Reads Mapped Confidently to Intronic Regions','Reads Mapped Confidently to Exonic Regions','Reads Mapped Antisense to Gene'), 
                        value = prettyNum(c(read_mapped_genome, read_uniq_mapped_genome, intergenic, intronic,exonic,antisense), big.mark = ','))
write.table(map_stats, paste0(opt$sample_id,"_map_stats.csv"),quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")

# ------------------------------ Calculate cell statistics

# calculation metrics about the barcoding and sequencing
gene_counts <- apply(filtered_mtx, 2, function(x) sum(x >= 1))
number_of_cell <- length(gene_counts)
mean_genes_cell <- round(mean(gene_counts))
med_genes_cell <- round(median(gene_counts))

# calculate mitochondrial percentage of genes
mito_genes <- grep("^MT-|mm10_mt-|GRCh38_MT-", rownames(filtered_mtx), value = TRUE)
percent_mito <- paste0(round(median(apply(filtered_mtx[mito_genes,names(gene_counts) ], 2, sum) / apply(filtered_mtx[,names(gene_counts)], 2, sum)),2),"%")


# Calculate the total counts for each cell
# Create a logical vector that indicates which genes have counts greater than 1
genes_with_counts_greater_than_1 <- rowSums(filtered_mtx) >= 1
# Calculate the total count for each cell
total_counts_per_cell <- rowSums(filtered_mtx[genes_with_counts_greater_than_1,names(gene_counts)])
# Calculate the mean counts per cell
mean_counts_per_cell <- round(mean(total_counts_per_cell))
# Calculate the median counts per cell
median_counts_per_cell <- round(median(total_counts_per_cell))


mean_read_cell <- round((total_number_reads/number_of_cell))
tot_genes_detected <- sum(rowSums(filtered_mtx)>=1)
cell_stats <- data.frame(stat = c('Estimated Number of Cells', 'Mean Reads per Cell', 'Mean Genes per Cell',
                                  'Median Genes per Cell', 'Total Genes Detected','Mean Counts per Cell',
                                  'Median Counts per Cell','Mitochondrial Ratio'), 
                        value = prettyNum(c(number_of_cell, mean_read_cell, mean_genes_cell,
                                        med_genes_cell, tot_genes_detected, mean_counts_per_cell,
                                        median_counts_per_cell, percent_mito
                                        ), big.mark = ','))
                                        
# ------------------------------ Write cell statistics
write.table(cell_stats, paste0(opt$sample_id,"_cell_stats.csv"),quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
