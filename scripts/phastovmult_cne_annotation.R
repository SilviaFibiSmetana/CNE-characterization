#Initialize the script
rm(list=ls())
args <- commandArgs ( TRUE )
if ( length ( args ) != 2 ){
  stop ( "Error. Check usage!" )
}
input <- args[1]
outputdir <- args[2]

if( ! file.exists(input) ){
  stop( paste(input, "does not exist!" ) )
}
outputdir <- paste(outputdir, "/", sep="")
if( ! file.exists(outputdir) ){
  stop( paste(outputdir, "does not exist!" ) )
}

sessionInfo()

##load the required packages
library("dplyr")
library("ggplot2")
library("tidyr")

#prepare output
output <- paste(outputdir,"/", basename(input), "_annotated", sep="")

# read input
tab <- read.table(input, sep="\t", header=FALSE)
head(tab,5)
dim(tab)

#first filter for the ones not overlapping with a gene (V7!=.) and set their overlap percentage to 0, calculate for the rest
tab1 <- tab %>%
  #filter(V5!="transcript") %>%
  mutate(V11=ifelse(V7==".", 0, 1)) %>%
  #mutate(V12=ifelse(V7==".", 0, 1) #%>%
  mutate(annotation=(ifelse(V11==1,"intronic", "intergenic"))) #%>%
  #  select(V1,GeneID, FixedCountExons,annotation) %>%
  #  separate(V1,c("V1","V2","V3"))
head(tab1,5)
dim(tab1)

tab2 <- tab1 %>%
  dplyr::select(V1, V2, V3, "annotation")
head(tab)
dim(tab2)

intro <- tab2 %>%
  filter(annotation == "intronic")
head(intro)
dim(intro)

inter <- tab2 %>%
  filter(annotation == "intergenic")
head(inter)
dim(inter)

write.table(tab2, paste(outputdir,"/", "noncoding_introinter.txt", sep=""), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)  
write.table(intro, paste(outputdir,"/", "noncoding_intro.txt", sep=""), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)  
write.table(inter, paste(outputdir,"/", "noncoding_inter.txt", sep=""), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)  
