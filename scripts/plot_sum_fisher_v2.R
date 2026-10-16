args <- commandArgs ( TRUE )
if ( length ( args ) != 2 ){
	stop ( "Error. Check usage!" )
}
inputfile <- args[1]
outputdir <- args[2]

if( ! file.exists(inputfile) ){
  stop( paste(inputfile, "does not exist!" ) )
}
if( ! file.exists(outputdir) ){
  stop( paste(outputdir, "does not exist!" ) )
}

library("dplyr")
library("tidyr")
library("ggplot2")
library("ggrepel")

#prepare output
output <- paste(outputdir,"/", basename(inputfile), "", sep="")

###################################################
### READ FILE
###################################################

df <- read.table(inputfile, sep = "\t", header = F)
head(df)
dim(df)

colnames(df)[colnames(df) == 'V2'] <- 'matchedCNEs'
colnames(df)[colnames(df) == 'V3'] <- 'totalCNEs'
colnames(df)[colnames(df) == 'V4'] <- 'SCREENelements'
colnames(df)[colnames(df) == 'V5'] <- 'pvalue'
colnames(df)[colnames(df) == 'V6'] <- 'ratio'
head(df)
dim(df)

df1 <- df %>%
  tidyr::separate(V1, c("clade", "a", "b", "c", "SCREENfile", "test"), sep = "_") %>%
  #separate(bla1, c("clade", "a", "b", "c"), sep = "_") %>%
  select(-a, -b, -c)
head(df1)
dim(df1)

  p <- ggplot(df1, aes(x=SCREENelements, y=matchedCNEs, label=SCREENfile)) +
    geom_point(aes(size = ratio, colour = pvalue)) +
    facet_wrap(~factor(clade, levels=c("Homininae", "Hominoidea", "Catarrhini", "Primates", "Euarchontoglires", "Boreoeutheria", "Eutheria", "Theria", "Mammals")), nrow = 3) + 
    geom_text_repel(aes(size = 3)) +
    ylab("overlapped CNEs") +
    scale_color_gradient(low="blue", high="red")

pdf(file=paste(output, "_", "SumFisher_dotplot.pdf", sep=""), width=9, height=8)
print(p)
dev.off()


#plot ratios stacked bar
p <- ggplot(df1, aes(fill = SCREENfile, y = ratio, x = factor(clade, levels=c("Homininae", "Hominoidea", "Catarrhini", "Primates", "Euarchontoglires", "Boreoeutheria", "Eutheria", "Theria", "Mammals")))) + 
  geom_bar(position = "stack", stat = "identity") + 
  #this could be added to add an asterisk for the not significant ones
  #geom_text(aes(label = ifelse(pvalue < 0.05, "", "*")), position = position_stack(vjust = 0.3), size = 6, color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(fill = "SCREENfile") +
  xlab("Clade")

pdf(file=paste(output, "_", "SumFisher_barchart.pdf", sep=""), width=9, height=8)
print(p)
dev.off()

