install.packages("BiocManager")

install.packages("optparse")

bio_pkgs <- c("rjson","R2HTML","DropletUtils","ggplot2")

BiocManager::install(bio_pkgs)
