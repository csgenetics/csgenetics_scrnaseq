install.packages("BiocManager")
install.packages("optparse")
install.packages("knitr")
install.packages("rmarkdown")
install.packages("dplyr")


bio_pkgs <- c("rjson","R2HTML","DropletUtils","ggplot2")

BiocManager::install(bio_pkgs)
