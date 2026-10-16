#Initialize the script
rm(list=ls())
args <- commandArgs(TRUE)
if ( length ( args ) != 3 ){
  stop ( "Error. Check usage!" )
}
cnes <- args[1]
random <- args[2]
outputdir <- args[3]

if( ! file.exists(cnes) ){
  stop( paste(cnes, "does not exist!" ) )
}
if( ! file.exists(random) ){
  stop( paste(random, "does not exist!" ) )
}
outputdir <- paste(outputdir, "/", sep="")
if( ! file.exists(outputdir) ){
  stop( paste(outputdir, "does not exist!" ) )
}

##load the required packages
library("dplyr")
library("tidyr")
library("ggplot2")

#prepare output
output <- paste(outputdir,"/",basename(cnes), "_corr", sep="")

##import the data
cne_table <- read.table(cnes, header = TRUE, sep = "\t")
#head(cne_table)
print(dim(cne_table))

rand <- read.table(random, header = TRUE, sep = "\t")
#head(rand)
print(dim(rand))

randommean <- mean(rand$my_mean)
#randommean

fold_changes <- cne_table %>%
  mutate(ratio=log2(my_mean/randommean)) %>%
  rename(my_mean_old = my_mean)
#head(fold_changes)
print(dim(fold_changes))

print(summary(fold_changes$ratio))

fold_changes.fixed <- fold_changes %>%
  mutate(my_mean=ifelse(as.numeric(ratio)<0, 0, ratio)) %>%
  select(-my_mean_old, -ratio)
#head(fold_changes.fixed)
print(dim(fold_changes.fixed))

print(summary(fold_changes.fixed$my_mean))

#export table
write.table(fold_changes.fixed %>% mutate(across(where(is.numeric), ~ round(., 4))), file=output, sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write(paste("# Output:", output), file="")
