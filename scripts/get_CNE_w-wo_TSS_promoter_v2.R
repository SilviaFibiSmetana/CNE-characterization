args <- commandArgs(TRUE)
if ( length(args) != 3 ){
  stop ( "Error. Check usage!" );
}
clade <- args[1];
input <- args[2];
outputdir <- args[3];

if( ! file.exists(input) ){
  stop( paste(input, "does not exist!" ) );
}
if( ! file.exists(outputdir) ){
  stop( paste(outputdir, "does not exist!" ) );
}

##load the required packages
library("data.table")
library("dplyr")
library("tidyr")

# import the data
tab <- read.table(input, header = F, sep = "\t")
print(dim(tab))

tab1 <- tab %>%
  group_by(V1, V2, V3) %>%  
  summarize(ID_combined = paste(unique(V7), collapse = ",")) 
print(dim(tab1))

output_file <- paste(outputdir, clade, "_estimated_CNEs_small_prom_overlap", sep = "") 
write.table(tab1, file=output_file, sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write(paste("# Output:", output_file), file="")

