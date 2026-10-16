rm(list=ls())
args <- commandArgs(TRUE)
if ( length ( args ) != 3 ){
  stop ( "Error. Check usage!" )
}
input <- args[1]
group <- args[2]
outputdir <- args[3]

if( ! file.exists(input) ){
  stop( paste(input, "does not exist!" ) )
}
if( ! file.exists(outputdir) ){
  stop( paste(outputdir, "does not exist!" ) )
}

library("data.table")
library("dplyr")
library("tidyr")

#prepare output
output <- paste(outputdir,"/", basename(group), "_estimated_CNE_numbers.txt", sep="")

############################################################
###Read input
############################################################

tab <- read.table(input, sep = "\t", header = T)
head(tab)
dim(tab)

tab1 <- cbind(RowNames= rownames(tab), tab)
head(tab1)
dim(tab1)


names(tab1)[2] ="Var"
dim(tab1)
head(tab1)

############################################################
###Process data
############################################################

tab2 <- tab1 %>%
  filter(Var == "1") %>%
  dplyr::select(RowNames)

dim(tab2)
head(tab2)

write.table(tab2, file=output, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write(paste("# Output:", output), file="")

