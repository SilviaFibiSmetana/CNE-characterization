rm(list=ls())
args <- commandArgs(TRUE)
if ( length ( args ) != 3 ){
  stop ( "Error. Check usage!" )
}
table <- args[1]
cne_conservation_file <- args[2]
dir <- args[3]

if( ! file.exists(table) ){
  stop( paste(table, "does not exist!" ) )
}
if( ! file.exists(cne_conservation_file) ){
  stop( paste(cne_conservation_file, "does not exist!" ) )
}
dir <- paste(dir, "/", sep="")
if( ! file.exists(dir) ){
  stop( paste(dir, "does not exist!" ) )
}

library("dplyr")
library("tidyr")
library("ggplot2")
library("gridExtra")
library("cowplot")
library("RColorBrewer")
library("fs")
library("data.table")
library("stringr")
library("scales")
library("colorspace")

#prepare output
output <- paste(dir, basename(table), "_out.txt", sep="")

#read file
consensus_table <- read.table(table, header=TRUE, check.names=FALSE, sep="\t")
head(consensus_table)
dim(consensus_table)

print(cne_conservation_file)
cne_conservation <- read.table(cne_conservation_file, header=TRUE, check.names=FALSE, sep="\t"); 
head(cne_conservation)
print(dim(cne_conservation))

# Fix the table so that coordinates are chromosome:start-end!!!
cne_conservation <- cne_conservation %>%
  separate(CNE, into=c("aux1", "aux2", "aux3"), sep="_") %>%
  unite(coor, aux1, aux2, sep=":") %>% 
  unite(coor, coor, aux3, sep="-")
# Rearrange
long.cne_conservation <- cne_conservation %>% 
  pivot_longer(!coor, names_to="clade", values_to="conserved") %>% 
  filter(!is.na(conserved)) %>% 
  select(-conserved)
# Keep only Mammals and Primates
long.cne_conservation <- long.cne_conservation %>%
  filter(clade %in% c("Mammals", "Primates"))
print(head(long.cne_conservation))
print(dim(long.cne_conservation))

#------------------------
# Pre-Processing
#------------------------


consensus_table[consensus_table == FALSE] <- 0
consensus_table[consensus_table == TRUE] <- 1
head(consensus_table)
consensus_table <- cbind(coor = rownames(consensus_table), consensus_table)
rownames(consensus_table) <- NULL
head(consensus_table)

# Add conservation
consensus_table.long <- consensus_table %>% pivot_longer(-coor, values_to="activity", names_to="tissue")
consensus_table.long <- consensus_table.long %>% left_join(long.cne_conservation, by="coor", relationship="many-to-many")
head(consensus_table.long)

#------------------------
# NUMBER OF TISSUES PER CNE
#------------------------
#Tissues Per CNEs
#Number of CNEs active in 0, 1, 2, etc tissues
NrOfActiveTissuesPerCNE <- consensus_table.long %>% group_by(clade, coor) %>% summarise(NrOfActiveTissues=sum(activity), NrOfTestedTissues=n_distinct(tissue)) %>% ungroup()

aux <- NrOfActiveTissuesPerCNE %>% group_by(clade, NrOfActiveTissues) %>% summarise(NrOfCNEs=n_distinct(coor)) %>% mutate(FractionOfCNEs=round(NrOfCNEs/sum(NrOfCNEs), digit=3))

head(aux)

#Number of CNEs active in at least one tissue
aux <- NrOfActiveTissuesPerCNE %>% group_by(clade) %>% summarise(NrOfActiveCNEs=sum(NrOfActiveTissues>0), FractionOfActiveCNEs=round(sum(NrOfActiveTissues>0)/n_distinct(coor), digit=3))
print(aux)

#Number of CNEs active in all tissues
aux <- NrOfActiveTissuesPerCNE %>% group_by(clade) %>% mutate(all=ifelse(NrOfActiveTissues==NrOfTestedTissues, TRUE, FALSE)) %>% summarise(NumberOfCNEsActiveInAllTissues=sum(all), FractionOfCNEsActiveInAllTissues=sum(all)/length(all))
print(aux)

#Mean and median number of tissues in which each CNE is active (Activity Level)
aux <- NrOfActiveTissuesPerCNE %>% group_by(clade) %>% summarise(MeanNrOfTissues=mean(NrOfActiveTissues), MedianNrOfTissues=median(NrOfActiveTissues))
print(aux)

#------------------------
# NUMBER OF CNES PER TISSUE
#------------------------

#Fraction of CNEs active in each tissue (Prevalence)
NrOfActiveCNEsPerTissue <- consensus_table.long %>% group_by(clade, tissue) %>% summarise(MeanNrOfActiveCNEs=mean(activity), MedianNumberOfActiveCNEs=median(activity), NrOfActiveCNEs=sum(activity), FractionOfActiveCNEs=round(sum(activity)/n_distinct(coor), digit=3)) %>% arrange(FractionOfActiveCNEs, .by_group=TRUE) %>% ungroup()
head(NrOfActiveCNEsPerTissue)

#------------------------
# Activity patterns
#------------------------
#HEatmap
# Reshape data
activity_matrix <- consensus_table.long %>% pivot_wider(names_from=tissue, values_from=activity)

clades <- consensus_table.long %>% select(clade) %>% pull() %>% unique()

write("# Number of activity patterns:", file="")
pattern_counts <- as.data.frame(activity_matrix) %>%
  group_by(across(-c(coor))) %>%
  summarise(count = n(), .groups = "drop") %>% 
  mutate(complexity = rowSums(select(., -c(clade, count))))

print(dim(pattern_counts))

# Exclude inactive CNEs
pattern_counts <- pattern_counts %>% filter(complexity != '0')
print(dim(pattern_counts))

#Theoretical model
p_avg <- consensus_table.long %>% group_by(clade) %>% summarise(MeanActivity=mean(activity), NrOfCNEs=n_distinct(coor))
print(p_avg)

nr_tissues <- consensus_table.long %>% summarise(n_distinct(tissue)) %>% pull()
print(nr_tissues)

pattern_counts <- pattern_counts %>% 
  left_join(p_avg, by="clade") %>% 
  mutate(ExpectedCount=MeanActivity**complexity * (1-MeanActivity)**(nr_tissues-complexity)*NrOfCNEs)


# Scatter plot: Number of tissues with activity in the pattern vs Number of elements with that pattern

p <- ggplot(pattern_counts, aes(x=complexity, y=count)) +
  # Layer the dots on top with jitter so they don't overlap
  geom_jitter(alpha=0.4, width=0.2, color="darkgray", size=2) + 
  stat_summary(fun=mean, geom="point",
               size=2, color = "red", shape = 18) + # Mean
  # Theoretical model
  stat_summary(fun=mean, geom="point", aes(y=ExpectedCount), 
               size=2, color = "black", shape = 22) + # Expected Value
  scale_y_log10(
    limits = c(0.1, NA),
    breaks = c(0.1, 1, 11, 101, 1001, 10001),
    labels = c("0", "1", expression(10^1), expression(10^2), expression(10^3), expression(10^4)),
    oob = scales::squish,
  ) +
  labs(x = "Pattern Complexity \n(Number of Tissues with Activity)", 
       y = "Number of CNEs") +
  #theme_clean() +
  theme(
    plot.title = element_text(hjust = 0.5),
  ) +
  facet_wrap(~clade)

output_pdf_file <- paste(output, "_Complexity.pdf", sep="")
pdf(output_pdf_file, width=12, height=5)
print(p)
dev.off()
