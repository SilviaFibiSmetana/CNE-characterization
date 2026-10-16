#Initialize the script
rm(list=ls())
args <- commandArgs(TRUE)
if ( length ( args ) != 5 ){
  stop ( "Error. Check usage!" )
}
prefix_cne_file <- args[1]
prefix_random_file <- args[2]
targetstring <- args[3]
cne_conservation_file <- args[4]
dir <- args[5]

prefix_cne_file <- "all_estimated_CNEs_"
prefix_random_file <- "random_estimated_CNEs_"
targetstring <- "_mean_corr" 
cne_conservation_file <- "output/02_mafoutput/summary_CNEs_in_clades.txt"
dir <- "output/03_functionaloutput/DNase"


if( ! file.exists(cne_conservation_file) ){
  stop( paste(cne_conservation_file, "does not exist!" ) )
}
dir <- paste(dir, "/", sep="")
if( ! file.exists(dir) ){
  stop( paste(dir, "does not exist!" ) )
}

library("dplyr")
library("tidyr")
library("fs")
library("stringr")
library("ggplot2")
library("scales")
library("purrr")

########################
# OUTPUT
########################
output <- paste(dir, prefix_cne_file, sep="")

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

# CNEs
pattern <- paste0(prefix_cne_file, ".*", targetstring, "$", sep = "") 
#pattern
myfiles <- fs::dir_ls(path = dir, regexp = pattern) 
if(length(myfiles)<=0){
	stop(paste("No files matching", pattern) )
}

epigenetic_signal_cnes <- NULL
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

  epigenetic_signal_cnes <- rbind(epigenetic_signal_cnes, ReadInFile)
}

# Generate a matrix with all CNEs and signals for all tissues
processed_epigenetic_signal_cnes <- epigenetic_signal_cnes %>%
  rename(tissue=File) %>%
  # Clean up tissue names
  mutate(tissue=basename(tissue)) %>%
  mutate(tissue=str_remove(tissue, prefix_cne_file)) %>%
  mutate(tissue=str_remove(tissue, targetstring)) %>%
  # Re-format coordinates
  unite("CNE", chr:start, sep = ":") %>%
  unite("CNE", CNE:end, sep = "-") %>%
  select(-C)

write("# CNEs:", file="")
print(head(processed_epigenetic_signal_cnes))
print(dim(processed_epigenetic_signal_cnes))

# Size-matched random genomic regions
pattern <- paste0(prefix_random_file, ".*", targetstring, "$", sep = "") 
#pattern
myfiles <- fs::dir_ls(path = dir, regexp = pattern) 
#myfiles_random 
if(length(myfiles)<=0){
	stop(paste("No files matching", pattern) )
}

epigenetic_signal_random <- NULL
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

  epigenetic_signal_random <- rbind(epigenetic_signal_random, ReadInFile)
}


# Generate a matrix with all random regions and signals for all tissues
processed_epigenetic_signal_random <- epigenetic_signal_random %>%
  rename(tissue=File) %>%
  # Clean up tissue names
  mutate(tissue=basename(tissue)) %>%
  mutate(tissue=str_remove(tissue, prefix_random_file)) %>%
  mutate(tissue=str_remove(tissue, targetstring)) %>%
  # Re-format coordinates
  unite("CNE", chr:start, sep = ":") %>%
  unite("CNE", CNE:end, sep = "-") %>%
  select(-C)

write("# Random genomic regions:", file="")
print(head(processed_epigenetic_signal_random))
print(dim(processed_epigenetic_signal_random))


#Data Pre-Processing
long.cne_signals <- processed_epigenetic_signal_cnes %>%
	rename(coor=CNE) %>%
	rename(Mean=my_mean) %>%
	mutate(type='foreground')
print(head(long.cne_signals))
print(dim(long.cne_signals))

long.random_signals <- processed_epigenetic_signal_random %>%
	rename(coor=CNE) %>%
	rename(Mean=my_mean) %>%
	mutate(type='background')

all_signals <- bind_rows(long.cne_signals, long.random_signals)
print(head(all_signals))
print(dim(all_signals))


# Calculating maximum signal across tissues for each CNE...
OneSignalPerCNE <- all_signals %>% 
	group_by(coor, type) %>% 
	slice_max(Mean, n=1, with_ties=FALSE) %>%
	arrange(desc(Mean)) %>%
	ungroup()
print(dim(OneSignalPerCNE))

# Adding Conservation Level
OneSignalPerCNE_with_conservation <- OneSignalPerCNE %>% 
	left_join(long.cne_conservation %>% 
		  rename(coor=CNE) %>%
		  mutate(type="foreground"), relationship="one-to-many", by=c("coor", "type"))
print(head(OneSignalPerCNE_with_conservation))

## Make sure the background has no conservation level!
print(OneSignalPerCNE_with_conservation %>% group_by(type, clade) %>% summarise(n_distinct(coor)))

# Enrichment Calculation
# Define thresholds
foreground_signal <- OneSignalPerCNE_with_conservation %>% 
  filter(type=="foreground") %>% 
  select(-clade) %>%
  unique() %>%
  select(Mean) %>% 
  pull()
print(length(foreground_signal))

intervals <- seq(0, 1, by=0.01)
thresholds <- quantile(foreground_signal, p=intervals)

# Map over thresholds
enrichment_results <- map_df(thresholds, function(t) {
  
  # 1. Calculate the global background proportion once for this threshold
  bg <- OneSignalPerCNE_with_conservation %>%
    filter(type == "background") %>%
    summarise(p = sum(Mean >= t), total=n())
  nr_bg_above <- bg %>% select(p) %>% pull()
  nr_bg_total <- bg %>% select(total) %>% pull()
  bg_prop <- nr_bg_above/nr_bg_total
  
  # 2. Calculate foreground stats per clade and divide by the global bg_prop
  OneSignalPerCNE_with_conservation %>%
    filter(type == "foreground") %>%
    group_by(clade) %>%
    summarise(
      threshold = t,
      NrForegroundAboveThreshold = sum(Mean >= t),
      TotalNrForeground = n(),
      # Proportion of elements above the threshold
      foreground_prop = sum(Mean >= t) / n(),
      NrBackgroundAboveThreshold = nr_bg_above,
      TotalNrBackground = nr_bg_total,
      bg_prop = bg_prop,
      # Use the pre-calculated bg_prop
      fold_enrichment = foreground_prop / bg_prop,
      SE = sqrt(1/NrForegroundAboveThreshold+1/NrBackgroundAboveThreshold),
      .groups = 'drop'
    )
})

#print(data.frame(enrichment_results))

# Figure
# Define the colors explicitly
col_primates <- "darkgreen"
col_mammals  <- "darkviolet"

# Prepare the data for plotting
plot_data <- enrichment_results %>%
  mutate(
    upper = 2**(log2(fold_enrichment) + (1.96 * SE)),
    lower = 2**(log2(fold_enrichment) - (1.96 * SE))
  )

## Filter out points where SE is too high (e.g., > 0.15)
plot_data <- plot_data %>%
  filter(SE < 0.15)

# Create the plot
p <- ggplot(plot_data, aes(x=foreground_prop, y=fold_enrichment, group=clade)) +
  # Area under the curve
  geom_area(alpha=0.2, position="identity", aes(fill=clade)) +
  # Add the main enrichment line
  geom_line(linewidth = 1, aes(color=clade)) +
  # Add points for clarity
  geom_point(size=1.5, aes(color=clade)) +
  # A dashed red line at 1 (no enrichment)
  geom_hline(yintercept=1, linetype="dashed", color="firebrick", alpha=0.5) +
  # Use a log10 scale for the y-axis
  scale_y_log10(
    limits = function(x) c(min(x, na.rm = TRUE), 10^ceiling(log10(max(x, na.rm = TRUE)))),
    breaks = function(limits) {
    # Automatically find integer powers of 10 that span the dynamic limits
    10^(floor(log10(limits[1])):ceiling(log10(limits[2])))
  },
  labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  #scale_y_log10(
  #  breaks = trans_breaks("log10", function(x) 10^x),
  #  labels = trans_format("log10", math_format(10^.x))
  #) +
  scale_x_continuous(
    labels = function(x) paste0(round((1 - x) * 100, 1), "%"),
    #            labels=scales::percent, 
  # Ensure the x-axis goes from 100% to 0% (lower threshold to higher threshols)
    transform="reverse", 
    #limits=c(0, NA), 
  ) +
  scale_fill_manual(values = c("Primates"=col_primates, "Mammals" = col_mammals), name = NULL) +
  scale_color_manual(values = c("Primates"=col_primates, "Mammals" = col_mammals), name = NULL) +
  # Formatting
  labs(
    x = "Signal Percentile (low to high intensity)",
    y = "Fold Enrichment",
  ) +
  theme(
    legend.position = c(0.15, 0.85),
    # Remove the background box and border around the legend
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank()
  )

output_pdf_file <- paste(output, "_enrichment.pdf", sep="")
pdf(output_pdf_file, width=5, height=5)
print(p)
dev.off()
