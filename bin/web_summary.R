#! /usr/bin/env Rscript

# load libraries
library(Matrix, quietly=T)
library(DropletUtils, quietly=T)
library(ggplot2, quietly=T)
library(scales, quietly=T)
library(rjson, quietly=T)
library(R2HTML, quietly=T)
library(optparse)

print_HTML <- function(seq_stats, map_stats ,cell_stats, dir, sample_id, plot){
  system(paste0('base64 ', plot ,' > ', dir, '/cell_caller.txt'))
  b64_bc <- readChar(paste0(dir, '/cell_caller.txt'), file.info(paste0(dir, '/cell_caller.txt'))$size)
  target <- HTMLInitFile("." ,filename=paste0(sample_id, '_summary'))
  HTML('<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=sans-serif">', file=target)
  HTML("<div class='title'>", file=target)
  HTML.title('QC Report per Sample', HR=1, file = target)
  HTML("</div>", file = target)
  HTML.title(sample_id, HR=3, file=target)
  HTML("<div id='wrapper'>", file=target) #start
  HTML("<div class='boxed' id='left' align='center'>", file=target)
  
  HTML('<table style="width:100%">', file=target)
  HTML(paste('<tr> <td align="center" style="font-size: 20px;" >', "Estimated Number of Cells", '</td> </tr> <tr> <td align="center" style="font-size: 18px;">', cell_stats$value[1] , '</td> </tr>'), file=target)
  HTML('</table>', file=target)
  
  HTML('<table style="width:100%">', file=target)
  HTML(paste('<tr> <td align="center" style="font-size: 20px;">', "Mean Reads per Cell", '</td> <td align="center" style="font-size: 20px;">', "Median Genes per Cell", '</td "> </tr>',
             '<tr> <td align="center" style="font-size: 20px;">', cell_stats$value[2], '</td> <td align="center" style="font-size: 20px;">', cell_stats$value[3], '</td> </tr>'), file=target)
  HTML('</table>', file=target)
  
  HTML.title('Sequencing', HR=2, file=target)
  HTML('<table style="width:100%">', file=target)
  HTML(paste('<tr> <td>', seq_stats$stat, '</td> <td align="right">', seq_stats$value, '</td>  </tr>'), file=target)
  HTML('</table>', file=target)
  
  HTML.title('Mapping', HR=2, file=target)
  HTML('<table style="width:100%">', file=target)
  HTML(paste('<tr> <td  >', map_stats$stat, '</td> <td  align="right" >', map_stats$value, '</td>  </tr>'), file=target)
  HTML('</table>', file=target)
  HTML("</div>", file = target)
  HTML("<div class='boxed' id='right' align='center'>", file=target)
  
  HTML('<table style="width:100%">', file=target)
  HTML(paste0("<img src='data:image/png;base64,", b64_bc, "' width=90%>"), file=target)
  HTML.title('Cells', HR=2, file=target)
  HTML(paste('<tr> <td>', cell_stats$stat, '</td> <td align="right">', cell_stats$value, '</td> </tr>'), file=target)
  HTML('</table>', file=target)

  HTML("</div>", file = target)
  HTML("</div>", file = target)
  HTML('<style type="text/css">
		.title {
    			background-color: #000000;
    			padding: 8px;
    			color: white;
    			top: 0;
    			left: 0;
    			z-index: 999;
    			width: 100%;
		}
		.boxed {
  			border: 1px solid #868D96;
  			padding: 10px;
  			margin: 10px;
		}
		h1 {
			font-family: "sans-serif";
			font-size: 33px;
                        margin: 3px;
		}
		h2 {
			font-family: "sans-serif";
			font-size: 26px;
			margin: 3px;
		}
		h3 {
			font-family: "sans-serif";
			font-size: 22px;
                        margin: 3px;
		}
		#wrapper {
  			display: flex;
		}
		#left {
  			width: 50%;
		}
		#right {
  			width: 50%;
		}
		table {
  			font-family: "sans-serif";
			font-size: 16px;
			border: 1px solid #868D96;
		}
		#mathplayer{
  			height: 0px;
		}
		</style>', file=target)
  HTMLEndFile()
}


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
    c('--plot'), action = 'store', type = 'character',
    help = 'Path to the cell caller plot  file. REQUIRED'),
  make_option(
    c('--cell_caller'), action = 'store', type = 'integer',
    help = 'The minimum number of nuclear genes detected for a barcode to be considered a cell. REQUIRED')
  )

opt <- parse_args(OptionParser(option_list=arg_list))

# Read in `matrix.mtx`
raw_mtx  <- readMM(opt$matrix)

# Read in `genes.tsv`
genes <- read.csv(opt$features, sep = '\t', header = F)
rownames(raw_mtx) <- genes[,1] # attach gene_ids

# Read in `barcodes.tsv`
cell_ids <- read.csv(opt$barcodes,sep = '\t', header = F)
colnames(raw_mtx) <- cell_ids$V1

out <- emptyDrops(raw_mtx) # get probability that each barcode is a cell
keep <- out$FDR <= 0.05 # define threshold probability for calling a cell
keep[is.na(keep)] <- FALSE
filt_mtx <- raw_mtx[,keep] # subset raw mtx to remove empty drops

write10xCounts(paste0('filter_count_matrix/', opt$sample_id), gene.symbol = genes[,2], filt_mtx, version = "3",overwrite=T) # write out filtered results



# load filter mtx
filt_mtx <- readMM(paste0('filter_count_matrix/', opt$sample_id,'/matrix.mtx.gz')) # load filtered mtx
# Read in filter `genes.tsv`
genes <- read.csv(paste0('filter_count_matrix/', opt$sample_id,'/features.tsv.gz'), sep = '\t', header = F)
rownames(filt_mtx) <- genes[,2] # attach gene_ids

# Read in filter `barcodes.tsv`
cell_ids <- read.csv(paste0('filter_count_matrix/', opt$sample_id,'/barcodes.tsv.gz'),sep = '\t', header = F)
colnames(filt_mtx) <- cell_ids$V1

### load multiqc json output 
data <- fromJSON(file=opt$multiqc_json)
total_number_reads   <- data$report_general_stats_data[[4]][[paste0(opt$sample_id,"_R1")]]$total_sequences

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


seq_stats <- data.frame(stat = c('Total Number of Reads', 'Valid Barcodes / CSGX Cell Barcodes', 'Sequencing Saturation', # get sequencing/alignment stats 
                                  'Q30 Bases in Barcode','Q30 Bases in RNA Read'), 
                        value = prettyNum(c(total_number_reads, valid_barcode, seq_sat, q30_barcode,q30_read), big.mark = ','))


map_stats <- data.frame(stat = c('Reads Mapped to Genome', 'Reads Mapped Confidently to Genome', 'Reads Mapped Confidently to Intergenic Regions',  
                                 'Reads Mapped Confidently to Intronic Regions','Reads Mapped Confidently to Exonic Regions','Reads Mapped Antisense to Gene'), 
                        value = prettyNum(c(read_mapped_genome, read_uniq_mapped_genome, intergenic, intronic,exonic,antisense), big.mark = ','))



#calculation metrics about the barcoding and sequencing process
ngenes <- opt$cell_caller
gene_counts <- apply(filt_mtx, 2, function(x) sum(x >= 1))
filter_genes_count <- gene_counts[gene_counts > ngenes]
number_of_cell <- length(filter_genes_count)
mean_genes_cell <- round(mean(filter_genes_count))
med_genes_cell <- round(median(filter_genes_count))


# calculate mitochondrial percentage of genes
mito_genes <- grep("^MT-|mm10_mt-|GRCh38_MT-", rownames(filt_mtx), value = TRUE)
percent_mito <- paste0(round(median(apply(filt_mtx[mito_genes,names(filter_genes_count) ], 2, sum) / apply(filt_mtx[,names(filter_genes_count)], 2, sum)),2),"%")


# Calculate the total counts for each cell
# Create a logical vector that indicates which genes have counts greater than 1
genes_with_counts_greater_than_1 <- rowSums(filt_mtx) >= 1
# Calculate the total count for each cell
total_counts_per_cell <- rowSums(filt_mtx[genes_with_counts_greater_than_1,names(filter_genes_count)])
# Calculate the mean counts per cell
mean_counts_per_cell <- round(mean(total_counts_per_cell))
# Calculate the median counts per cell
median_counts_per_cell <- round(median(total_counts_per_cell))


mean_read_cell <- round((total_number_reads/number_of_cell))
tot_genes_detected <- sum(rowSums(filt_mtx)>=1)
cell_stats <- data.frame(stat = c('Estimated Number of Cells', 'Mean Reads per Cell', 'Mean Genes per Cell',
                                  'Median Genes per Cell', 'Total Genes Detected','Mean Counts per Cell',
                                  'Median Counts per Cell','Mitochondrial Ratio'), 
                        value = prettyNum(c(number_of_cell, mean_read_cell,mean_genes_cell,
                                        med_genes_cell, tot_genes_detected,mean_counts_per_cell,median_counts_per_cell,percent_mito
                                        ), big.mark = ','))

combine_df <- rbind(seq_stats,map_stats,cell_stats)
colnames(combine_df) <- c("scRNA_Metrics",opt$sample_id)
# Save the final result
write.csv(combine_df, paste0(opt$sample_id,"_scRNA_Metrics.csv"), row.names = FALSE)

print_HTML(seq_stats = seq_stats, map_stats= map_stats ,cell_stats = cell_stats, dir = paste0('filter_count_matrix/', opt$sample_id) , sample_id = opt$sample_id, plot = opt$plot) # output a HTML summary of the run

print('Summary Report Done!')

