#Initialize the script
rm(list=ls())
args <- commandArgs(TRUE)
if ( length ( args ) != 2 ){
    stop ( "Error. Check usage!" )
}
input <- args[1]
outputdir <- args[2]; 

if( !file.exists(input) ){
  stop( paste(input, "does not exist!" ) )
}
outputdir <- paste(outputdir, "/", sep="")
if( !file.exists(outputdir) ){
  stop( paste(outputdir, "does not exist!" ) )
}

##load the required packages
library("dplyr")
library("tidyr")
library("ggplot2")

#prepare output
output <- paste(outputdir,basename(input), "_fccategory", sep="")

##import the data
tab <- read.table(input, header = TRUE, sep = "\t")
head(tab)
dim(tab)


tab1 <- tab %>%
  mutate(category = ifelse(my_mean < 0.5, "1",
                           ifelse(my_mean >= 0.5 & my_mean < 1, "2",
                                  ifelse(my_mean >= 1 & my_mean < 1.5 , "3", 
                                         ifelse(my_mean >= 1.5 & my_mean < 2 , "4", "5"))))) 

head(tab1)
dim(tab1)

tab2 <- tab1 %>%
  group_by(category) %>%
  summarize(n=n())
head(tab2)
dim(tab2)

#export table
write.table(tab1,file=output, sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write(paste("# Output", output), file="")
