rm(list=ls())
args <- commandArgs ( TRUE )
if ( length ( args ) != 4 ){
	stop ( "Error. Check usage!" )
}
alltab <- args[1]
mamtab <- args[2]
pritab <- args[3]
outputdir <- args[4]

if( ! file.exists(alltab) ){
  stop( paste(alltab, "does not exist!" ) )
}
if( ! file.exists(mamtab) ){
  stop( paste(mamtab, "does not exist!" ) )
}
if( ! file.exists(pritab) ){
  stop( paste(pritab, "does not exist!" ) )
}
if( ! file.exists(outputdir) ){
  stop( paste(outputdir, "does not exist!" ) )
}

library("dplyr")
library("tidyr")

#prepare output
output <- paste(outputdir,"/", "summary_CNEs_in_clades.txt", sep="")

############################################################
###Read input
############################################################

all <- read.table(alltab, sep = "\t", header = F)
head(all)
dim(all)

mam <- read.table(mamtab, sep = "\t", header = F)
head(mam)
dim(mam)

pri <- read.table(pritab, sep = "\t", header = F)
head(pri)
dim(pri)

############################################################
###Combine
############################################################

all1 <- all %>%
  unite("cne", V1:V3)
names(all1)[names(all1) == "cne"] <- "V1"
head(all1)
dim(all1)

tab <- all1 %>% left_join(mam %>% mutate(Mammals = 1), by = "V1")
head(tab)
dim(tab)

tab1 <- tab %>% left_join(pri %>% mutate(Primates = 1), by = "V1")
names(tab1)[names(tab1) == "V1"] <- "CNE"
head(tab1)
dim(tab1)

write.table(tab1, file=output, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write(paste("# Output:", output), file="")
