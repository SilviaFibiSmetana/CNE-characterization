#Initialize the script
rm(list=ls())
args <- commandArgs(TRUE)
if ( length ( args ) != 4 ){
  stop ( "Error. Check usage!" )
}
clade <- args[1]
targetstring <- args[2]
cne_conservation_file <- args[3]
dir <- args[4]

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
output <- paste(dir, clade, sep="")

########################
# INPUT
########################

# Read CNEs with the clades in which they are conserved
cne_conservation <- read.table(cne_conservation_file, header=TRUE, check.names=FALSE, sep="\t") 
# Fix CNE format
cne_conservation <- cne_conservation %>%
	mutate(CNE=sub("^([^_]+)_([^_]+)_(.*)$", "\\1:\\2-\\3", CNE))
print(dim(cne_conservation))

# Rearrange
long.cne_conservation <- cne_conservation %>% 
  pivot_longer(!CNE, names_to="clade", values_to="conserved") %>% 
  filter(!is.na(conserved)) %>% 
  select(-conserved)

# Keep only Mammals and Primates
long.cne_conservation <- long.cne_conservation %>%
  filter(clade %in% c("Mammals", "Primates"))
print(head(long.cne_conservation))


# Import files with DNase/H3K27ac signals
pattern <- paste0(clade, ".*", targetstring, "$", sep = "") 
#pattern
myfiles <- fs::dir_ls(path = dir, regexp = pattern) 
#myfiles 

epigenetic_signal_table <- NULL
for(filename in myfiles){ 
  if( !file.exists(filename) ){
    stop( paste(filename, "does not exist!" ) )
  }

  ReadInFile <- read.table(file=filename, header=T, na.strings="NULL")
  ReadInFile$File <- filename
  
  aux <- ReadInFile$File
  parts <- strsplit(aux, "_")[[1]]
  aux <- tail(parts, 2)
  aux <- paste(aux, collapse = "_")
  ReadInFile$C <- aux

  epigenetic_signal_table <- rbind(epigenetic_signal_table, ReadInFile)
}

# Generate a matrix with all CNEs and signals for all tissues
processed_epigenetic_signal_table <- epigenetic_signal_table %>%
  rename(tissue=File) %>%
  # Clean up tissue names
  mutate(tissue=basename(tissue)) %>%
  mutate(tissue=str_remove(tissue, clade)) %>%
  mutate(tissue=str_remove(tissue, targetstring)) %>%
  # Re-format coordinates
  unite("CNE", chr:start, sep = ":") %>%
  unite("CNE", CNE:end, sep = "-")

# Average signal across random regions corresponding to the same CNE
processed_epigenetic_signal_table <- processed_epigenetic_signal_table %>%
  group_by(tissue) %>%
  mutate(randommean=mean(my_mean)) %>% 
  ungroup() %>%
  mutate(ratio=log2(as.numeric(my_mean)/randommean))

processed_epigenetic_signal_table <- processed_epigenetic_signal_table %>%
	rename(my_old_mean=my_mean)

processed_epigenetic_signal_table <- processed_epigenetic_signal_table %>%
  mutate(my_mean=ifelse(as.numeric(ratio)<0, 0, ratio)) %>%
  select(-my_old_mean, -ratio)

processed_epigenetic_signal_table.wide <- processed_epigenetic_signal_table %>%
	select(CNE, tissue, my_mean) %>%
  pivot_wider(names_from = tissue, values_from = my_mean)

## WRITE TABLE TO A FILE
write.table(processed_epigenetic_signal_table.wide %>% mutate(across(where(is.numeric), ~ round(., 4))), quote=FALSE, sep="\t", row.names=FALSE, file=paste(output, "Signal.tsv", sep="_"))
write(paste("# Output", paste(output, ".tsv", sep="")), file="")

n_cnes <- nrow(processed_epigenetic_signal_table.wide)
n_tissues <- processed_epigenetic_signal_table.wide %>% select(-CNE) %>% ncol()


########################
# COMPUTE SIGNAL CATEGORIES
########################

long.cne_signals <- processed_epigenetic_signal_table.wide %>% 
	pivot_longer(!CNE, names_to="tissue", values_to="LogRatio")

# Compute summary statistics
avg.cne_signals <- long.cne_signals %>% 
	group_by(tissue) %>% 
	summarise(mean(LogRatio), median(LogRatio), min(LogRatio), max(LogRatio))
print(head(avg.cne_signals))

# Compute signal categories
long.cne_signals <- long.cne_signals %>% 
  mutate(Category=ifelse(LogRatio<0.5, 1, 
                         ifelse(LogRatio<1, 2, 
                                ifelse(LogRatio<1.5, 3, 
                                       ifelse(LogRatio<2, 4, 5)))))
n_categories <- long.cne_signals %>% 
	select(Category) %>% 
	unique() %>% 
	nrow()

########################
# ADD CONSERVATION CLADE
########################

long.cne_signals_with_conservation <- long.cne_signals %>% 
  left_join(long.cne_conservation, relationship="many-to-many", by="CNE")
head(long.cne_signals_with_conservation)

#Total Number of CNEs
aux <- long.cne_signals_with_conservation %>% group_by(clade) %>% summarise(NrCNEs=n_distinct(CNE))
head(data.frame(aux))

########################
# Number of CNEs per Category, per Tissue
########################

## Pre-calculate fractions
df_plot <- long.cne_signals_with_conservation %>%
  group_by(clade, tissue, Category) %>%
  summarise(NrOfCNEs=n_distinct(CNE)) %>% 
  group_by(clade, tissue) %>%
  mutate(FractionOfCNEs = round(NrOfCNEs/sum(NrOfCNEs), digit=3)) %>%
  mutate(clade=substr(clade, 0, 1))
head(data.frame(df_plot))

p <- ggplot(df_plot, aes(x=clade, y=FractionOfCNEs, fill=factor(Category))) + 
  geom_col() +
  labs(title='', x='', y='Fraction of CNEs') +
  theme_bw() +
  theme(legend.position="top",
        text = element_text(color="black"),
        axis.text = element_text(color="black"),
        panel.spacing = unit(0.1, "lines"),
        strip.text = element_text(size=6, face="bold")
  ) +
  scale_fill_discrete_sequential(palette = "Reds") +
  guides(fill = guide_legend(title = "Category")) +
  facet_wrap(~tissue, ncol=n_tissues)

output_pdf_file <- paste(output, "PercentInCategories.pdf", sep="_")
pdf(output_pdf_file, width=12, height=5)
print(p)
dev.off()

########################
# NUMBER OF ACCESSIBLE/ACETYLATED CNES
########################

active_categories <- c(4,5)
long.cne_signals_with_conservation <- long.cne_signals_with_conservation %>% 
  mutate(active=ifelse(Category %in% active_categories, TRUE, FALSE))
head(long.cne_signals_with_conservation)

#Total number of accessible/acetylated CNEs
NrOfTissuesShowingActivity <- long.cne_signals_with_conservation %>% 
  group_by(clade, CNE) %>% 
  summarise(NrActiveCNEs=sum(active)) %>%
  ungroup()

df_plot <- NrOfTissuesShowingActivity

# Define the colors explicitly
col_primates <- "darkgreen"
col_mammals  <- "darkviolet"
p <- ggplot(df_plot, aes(x=NrActiveCNEs, group=clade, fill=clade)) +
  geom_histogram(aes(y = after_stat(density)),
                 color=NA, bins=19, position="dodge", alpha=0.5)+
  scale_fill_manual(values = c("Primates"=col_primates, "Mammals" = col_mammals), name=NULL) +
  labs(title="",x="Number of Tissues",  y="Relative Frequency")+
  theme_classic()

output_pdf_file <- paste(output, "NrOfTissuesShowingActivity.pdf", sep="_")
pdf(output_pdf_file, width=8, height=8)
print(p)
dev.off()

#CNEs that are accessible/acetylated in at least one tissue
aux <- NrOfTissuesShowingActivity %>% mutate(active=NrActiveCNEs>0) %>% group_by(clade, active) %>% summarise(NrCNEs=n_distinct(CNE)) %>% group_by(clade) %>% mutate(FractionOfCNEs=round(NrCNEs/sum(NrCNEs), digit=3))
head(data.frame(aux))

#CNEs that are accessible/acetylated in exactly one tissue
aux <- NrOfTissuesShowingActivity %>% mutate(aux=(NrActiveCNEs==1)) %>% group_by(clade) %>% summarise(NrOfCNEs=sum(aux), FractionOfCNEs=round(sum(aux)/n_distinct(CNE), digit=3))
head(data.frame(aux))

#CNEs that are accessible/acetylated in three or more tissues
aux <- NrOfTissuesShowingActivity %>% mutate(aux=(NrActiveCNEs>=3)) %>% group_by(clade) %>% summarise(NrOfCNEs=sum(aux), FractionOfCNEs=round(sum(aux)/n_distinct(CNE), digit=3))
head(data.frame(aux))

#CNEs that are accessible/acetylated in all tissues
aux <- NrOfTissuesShowingActivity %>% mutate(aux=(NrActiveCNEs==n_tissues)) %>% group_by(clade) %>% summarise(NrOfCNEs=sum(aux), FractionOfCNEs=round(sum(aux)/n_distinct(CNE), digit=3))
head(data.frame(aux))

#Median and mean number of tissues in which each CNE is accessible/acetylated (including non-accessible/acetylated CNEs)
print(NrOfTissuesShowingActivity %>% group_by(clade) %>% summarise(median(NrActiveCNEs), mean(NrActiveCNEs)))

#Median and mean number of tissues in which each accessible/acetylated CNE is accessible/acetylated
print(NrOfTissuesShowingActivity %>% group_by(clade) %>% filter(NrActiveCNEs>0) %>% summarise(median(NrActiveCNEs), mean(NrActiveCNEs)))

#Number of accessible/acetylated CNEs, per tissue
ActiveCNEs <- long.cne_signals_with_conservation %>% group_by(clade, tissue) %>% 
  summarise(NrActiveCNEs=sum(active), FractionActiveCNEs=round(sum(active)/n_distinct(CNE), digit=3)) %>% 
  group_by(clade) %>%
  arrange(desc(FractionActiveCNEs), .by_group=TRUE)
head(data.frame(ActiveCNEs))

## REARRANGE AND WRITE TABLE TO A FILE
all_signals.wide <- long.cne_signals %>% 
	select(CNE, tissue, Category) %>% 
	pivot_wider(names_from=tissue, values_from=Category)
write.table(all_signals.wide, quote=FALSE, sep="\t", row.names=FALSE, file=paste(output, "Categories.tsv", sep="_"))
write(paste("# Output", paste(output, ".tsv", sep="")), file="")

