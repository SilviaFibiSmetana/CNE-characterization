rm(list=ls())
args <- commandArgs ( TRUE )
if ( length ( args ) != 4 ){
  stop ( "Error. Check usage!" )
}
alltab <- args[1]
cladetab <- args[2]
clade <- args[3]
outputdir <- args[4]

if( ! file.exists(alltab) ){
  stop( paste(alltab, "does not exist!" ) )
}
if( ! file.exists(cladetab) ){
  stop( paste(cladetab, "does not exist!" ) )
}
outputdir <- paste(outputdir, "/", sep="")
if( ! file.exists(outputdir) ){
  stop( paste(outputdir, "does not exist!" ) )
}

library("dplyr")
library("tidyr")

#prepare output
x <- basename(alltab)
x
y <- sub("^all", "", x)
y
output <- paste(outputdir,"/", clade, y, sep="")

############################################################
###Read input
############################################################

all <- read.table(alltab, sep = "\t", header = T)
head(all)
dim(all)

tab <- read.table(cladetab, sep = "\t", header = F)
head(tab)
dim(tab)

############################################################
###Combine
############################################################

all1 <- all %>%
  unite("CNE", chr:end, sep = "_")
head(all1)
dim(all1)

tab1 <- tab %>%
  unite("CNE", V1:V3, sep = "_")
head(tab1)
dim(tab1)

spec <- all1 %>% inner_join(tab1, by = "CNE")
head(spec)
dim(spec)

spec1 <- spec %>%
  separate("CNE", c("chr", "start", "end"))
head(spec1)
dim(spec1)

write.table(spec1, file=output, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write(paste("# Output", output), file="")
