rm(list=ls())
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

if (ncol(df) != 6) {
  stop("Error: Expected 6 columns, but found ", ncol(df), ".")
}

# Rename columns
df <- df %>%
  rename(
    label          = V1,
    matchedCNEs    = V2,
    totalCNEs      = V3,
    SCREENelements = V4,
    pvalue         = V5,
    ratio          = V6
)

df <- df %>% 
	mutate(label=gsub("_estimated_CNEs_AND", "", label))
print(head(df))

df1 <- df %>% 
	separate_wider_delim(label, names=c("clade", "Tissueshort", "element", "test"), delim="_")
print(head(df1))

# Filter elements
df1.filtered <- df1 %>% 
	filter(element %in% c("pELS", "dELS", "PLS"))

# Visualize
p <- ggplot(df1.filtered, aes(x=SCREENelements, y=matchedCNEs, label=Tissueshort, group=clade)) +
      geom_point(aes(size=ratio, colour=clade)) +
      facet_wrap(~element, scales = "free") + 
      geom_text_repel(aes(size = 3)) +
      ylab("Number of CNEs") +
      xlab("Number of SCREEN elements") 
  
pdf(file=paste(output, "_dotplot.pdf", sep=""), width=10, height=5)
print(p)
dev.off();  

p <- ggplot(df1.filtered, aes(fill = Tissueshort, y = ratio, x = clade)) + 
  geom_bar(position = "stack", stat = "identity") + 
  facet_wrap(~element, scales = "free") +
  #this could be added to add an asterisk for the not significant ones
  #geom_text(aes(label = ifelse(pvalue < 0.05, "", "*")), position = position_stack(vjust = 0.3), size = 6, color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(fill = "Tissue") +
  xlab("Clade")

pdf(file=paste(output, "_barchart.pdf", sep=""), width=10, height=5)
print(p)
dev.off()

