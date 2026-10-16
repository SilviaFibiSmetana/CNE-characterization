#Initialize the script
rm(list=ls())
args <- commandArgs(TRUE)
if ( length ( args ) != 3 ){
    stop ( "Error. Check usage!" )
}
file1 <- args[1]
file2 <- args[2]
dir <- args[3]

if( ! file.exists(file1) ){
  stop( paste(file1, "does not exist!" ) )
}
if( ! file.exists(file2) ){
  stop( paste(file2, "does not exist!" ) )
}
dir <- paste(dir, "/", sep="")
if( ! file.exists(dir) ){
  stop( paste(dir, "does not exist!" ) )
}

library("dplyr")
library("tidyr")
library("ggplot2")

#prepare output
output <- paste(dir, basename(file1), "_consensus", sep="")
output

#import files
taba <- read.table(file1, header = T)
head(taba)
dim(taba)

tabb <- read.table(file2, header = T)
head(tabb)
dim(tabb)

#preporcessing
taba1 <- taba %>% 
  pivot_longer(!CNE, names_to = "tissues", values_to = "dnase_category") %>%
  unite("A", CNE:tissues)
head(taba1)
dim(taba1)

tabb1 <- tabb %>% 
  pivot_longer(!CNE, names_to = "tissues", values_to = "h3k27ac_category") %>%
  unite("A", CNE:tissues)
head(tabb1)
dim(tabb1)

tabcon <- inner_join(taba1, tabb1, by = "A")
head(tabcon)
dim(tabcon)

tabcon1 <- tabcon %>%
  mutate(consensus=ifelse((dnase_category==5 | dnase_category==4) & (h3k27ac_category==4 | h3k27ac_category==5), "both", 
                          ifelse((dnase_category==5 | dnase_category==4) & h3k27ac_category<4, "DNase", 
                                 ifelse(dnase_category<4 & (h3k27ac_category==4 | h3k27ac_category==5), "H3K27ac", "low"))))
head(tabcon1)
dim(tabcon1)

tabout <- tabcon1 %>%
  separate(A, c("CNE", "tissues"), sep = "_") %>%
  select(CNE, tissues, consensus)
head(tabout)
dim(tabout)

write.table(tabout, file = paste(output, "list", sep=""), sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write(paste("# Output:", paste(output, "list", sep="")), file="")

tabout1 <- tabout %>%
  pivot_wider(names_from = tissues, values_from = consensus)
tabout1 <- as.data.frame(tabout1)
rownames(tabout1) <- tabout1[,1]
tabout1$CNE <- NULL
tabout1[] <- lapply(tabout1, function(x) ifelse(x == "both", "TRUE", "FALSE"))
head(tabout1)
dim(tabout1)

write.table(tabout1, file = paste(output, "activelist", sep=""), sep="\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
write(paste("# Output:", paste(output, "activelist", sep="")), file="")


tabcon2 <- tabcon1 %>%
  group_by(consensus) %>%
  summarise(n=n())
head(tabcon2)
dim(tabcon2)

write.table(tabcon2, file = output, sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write(paste("# Output:", output), file="")

tabcondnase5 <- tabcon %>% 
  filter(dnase_category==4 | dnase_category==5) %>%
  group_by(h3k27ac_category) %>%
  summarise(n=n())
head(tabcondnase5)
dim(tabcondnase5)

tabconh3k27ac45 <- tabcon %>% 
  filter(h3k27ac_category==4 | h3k27ac_category==5) %>%
  group_by(dnase_category) %>%
  summarise(n=n())
head(tabconh3k27ac45)
dim(tabconh3k27ac45)

write.table(tabcondnase5, file = paste(output, "_dnase5.txt", sep=""), sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(tabconh3k27ac45, file = paste(output, "_h3k27ac45.txt", sep=""), sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
