#Initialize the script
rm(list=ls())
args <- commandArgs(TRUE)
if ( length ( args ) != 2 ){
  stop ( "Error. Check usage!" )
}
file <- args[1]
dir <- args[2]

if( ! file.exists(file) ){
  stop( paste(file, "does not exist!" ) )
}
dir <- paste(dir, "/", sep="")
if( ! file.exists(dir) ){
  stop( paste(dir, "does not exist!" ) )
}

library("dplyr")
library("tidyr")
library("ggplot2")

#prepare output
output <- paste(dir, basename(file), "_corr", sep="")
output

#import file
tab <- read.table(file, header = T)
tab <- rename(tab, my_mean_old = my_mean)
head(tab)
dim(tab)

tab1 <- tab %>%
  mutate(my_mean=ifelse(as.numeric(my_mean_old)<0, 0, my_mean_old)) %>%
  select(-my_mean_old)
head(tab1)
dim(tab1)

write.table(tab1, file = output, sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write(paste("# Output", output), file="")
