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

#first filter for the ones not overlapping with a CDS (V3==.) and set their overlap percentage to 0, calculate for the rest
tab1 <- tab %>%
  mutate(V8=ifelse(V4==".", 0, (V7/V2))) %>%
  mutate(V9=ifelse(V4==".", 0, (V2/V7))) 
head(tab1,5)
dim(tab1)


tab2 <- tab1 %>%
  group_by(V1) %>% 
  summarize(GeneID=paste(as.character(V6),collapse="/"), TotPerc=sum(V8), ExonIncl=sum(V9), CountExons=n()) %>%
  mutate(FixedCountExons=ifelse(GeneID==".", 0, CountExons)) %>%
  mutate(annotation=ifelse(TotPerc==0,"noncoding", ifelse(TotPerc>=0.3,"coding", ifelse((TotPerc<0.3) & (ExonIncl>=1), "coding", "ambiguous")))) %>%
  dplyr::select(V1,GeneID, FixedCountExons,annotation) %>%
  separate(V1,c("V1","V2","V3"))
head(tab2,5)
dim(tab2)


#write table
write.table(tab2, output, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)  
write(paste("# Output:", output), file="")
