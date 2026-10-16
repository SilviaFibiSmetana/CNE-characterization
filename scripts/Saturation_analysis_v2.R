#Initialize the script
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
output <- paste(dir, basename(table), "_saturation.txt", sep="")

#------------------------
# Parameters
#------------------------

set.seed(123) 

n_perms <- 100

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

n_tissues <- ncol(consensus_table)-1
n_tissues

## Add conservation
consensus_table_WithConservation <- consensus_table %>% left_join(long.cne_conservation, by="coor")
head(consensus_table_WithConservation)

#------------------------
# Functions
#------------------------

# Get one random accumulation path
get_path <- function(m) {
  m_shuffled <- m[, sample(ncol(m))]
  # Check if any column up to 'i' has a 1 (rowSums > 0)
  sapply(1:ncol(m_shuffled), function(i) {
    sum(rowSums(m_shuffled[, 1:i, drop = FALSE]) > 0)
  })
}

GetSaturationModel <- function(mat, n_perms){
  mat$coor <- NULL
  #print(head(mat))
  mat <- data.matrix(mat)
  print(dim(mat))
  
  n_cols <- ncol(mat)
  
  # Run permutations and calculate the mean/SD
  results <- replicate(n_perms, get_path(mat))
  
  plot_df <- data.frame(
    Columns = 1:n_cols,
    Median_Richness = apply(results, 1, median),
    Mean_Richness = rowMeans(results),
    SD_Richness = apply(results, 1, sd)
  )
  
  # Fit a Model
  ## Fit a Linear model on reciprocals
  lin_fit <- lm(1/Mean_Richness ~ I(1/Columns), data = plot_df)
  ## Extract coefficients
  intercept <- coef(lin_fit)[1]
  slope <- coef(lin_fit)[2]
  ## Calculate better starting guesses
  start_Vmax <- 1 / intercept
  start_Km <- slope * start_Vmax
  ## Perform a non-linear least squares regression to fit the data to a Michaelis-Menten model:
  fit <- nls(Mean_Richness ~ (Vmax * Columns) / (Km + Columns),
             data = plot_df,
             start = list(Vmax=start_Vmax, Km=start_Km))
  print(coef(fit))
  
  # Add Model to the data frame
  plot_df$ModelPrediction <- predict(fit, plot_df)
  
  # Analyze slope (saturation)
  v_max <- coef(fit)["Vmax"]
  k_m <- coef(fit)["Km"]
  N <- ncol(mat)
  final_slope <- as.numeric((v_max * k_m) / (k_m + N)^2)
  
  saturation <- list(fit=fit, plot_df=plot_df, final_slope=final_slope)
  
  return(saturation)
}


VisualizeSaturation <- function(plot_df){
  
  p <- ggplot(plot_df, aes(x=Columns, y=Mean_Richness)) +
    # Uncertainty ribbon
    geom_ribbon(aes(ymin = Mean_Richness - SD_Richness,
                    ymax = Mean_Richness + SD_Richness),
                fill = "steelblue", alpha = 0.2) +
    # Raw average line
    geom_line(color="steelblue", linewidth = 1) +
    geom_line(aes(y=ModelPrediction), color="firebrick", linetype="dashed") +
    labs(x = "Number of Tissues Sampled",
         y = "Cumulative Number of\n Active CNEs") +
    theme(
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 14), 
      strip.text = element_text(color = "black", size = 12),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
    ) +
    facet_wrap(~clade, scales="free_y")
  
  return(p)
}

#------------------------
# Saturation Analysis
#------------------------

saturation_result <- consensus_table_WithConservation %>%
  group_by(clade) %>%
  # Wrap the function output in a list so reframe can handle it
  reframe(res = list(GetSaturationModel(pick(everything()), n_perms)))

all_plot_df <- list()
for (i in 1:nrow(saturation_result)) {
  tmp_df <- saturation_result$res[[i]][['plot_df']]
  tmp_df$clade <- saturation_result$clade[i]
  all_plot_df[[i]] <- tmp_df
}
all_plot_df <- bind_rows(all_plot_df)

annotation <- saturation_result %>% 
  rowwise() %>% 
  mutate(Vmax=round(coef(res[["fit"]])["Vmax"], digit=0), Km=coef(res[["fit"]])["Km"]) %>%
  mutate(model = paste0("y == frac(", round(Vmax, digit=1), " %.% x, ", round(Km, digit=1), " + x)"))
annotation$x_pos = -Inf
annotation$y_pos = Inf

aux <- consensus_table_WithConservation %>% group_by(clade) %>% summarise(NrCNEs=n_distinct(coor))
annotation <- annotation %>% left_join(aux, by="clade")

p <- VisualizeSaturation(all_plot_df)

p <- p +
  geom_text(
    data = annotation,
    aes(x = x_pos, y = y_pos, label = model),
    parse = TRUE,        
    inherit.aes = FALSE, 
    hjust = -0.2,    
    vjust = 1.2,    
    size = 3,
    family = "sans"
  )

output_pdf_file <- paste(output, ".pdf", sep="")
pdf(output_pdf_file, width=12, height=5)
print(p)
dev.off()

#Mammalian subset
results <- saturation_result %>% filter(clade == "Mammals") %>% select(res)
fit <- results$res[[1]]$fit
final_slope <- results$res[[1]]$final_slope
saturation_pct <- (results$res[[1]]$plot_df$Mean_Richness[n_tissues] / coef(fit)["Vmax"]) * 100

NrCNEs <- annotation %>% filter(clade=="Mammals") %>% select(NrCNEs) %>% pull()

# Calculate the slope at the current sample size
print(final_slope)

print(saturation_pct)

n_active_cnes <- as.numeric(coef(fit)["Vmax"]*n_tissues/(coef(fit)["Km"]+n_tissues))
print(round(n_active_cnes, digit=0))

print(round(n_active_cnes/NrCNEs, digit=4))

#(Mean) Number of active CNEs at saturation
n_active_cnes_at_saturation <- as.numeric(coef(fit)["Vmax"])
print(round(n_active_cnes_at_saturation, digit=0))

print(round(n_active_cnes_at_saturation/NrCNEs, digit=4))

#Primate subset
results <- saturation_result %>% filter(clade == "Primates") %>% select(res)
fit <- results$res[[1]]$fit
final_slope <- results$res[[1]]$final_slope
saturation_pct <- (results$res[[1]]$plot_df$Mean_Richness[n_tissues] / coef(fit)["Vmax"]) * 100

NrCNEs <- annotation %>% filter(clade=="Primates") %>% select(NrCNEs) %>% pull()

# Calculate the slope at the current sample size
print(final_slope)

print(saturation_pct)

#(Mean) Number of active CNEs for 19 tissues
n_active_cnes <- as.numeric(coef(fit)["Vmax"]*n_tissues/(coef(fit)["Km"]+n_tissues))
print(round(n_active_cnes, digit=0))

print(round(n_active_cnes/NrCNEs, digit=4))

n_active_cnes_at_saturation <- as.numeric(coef(fit)["Vmax"])
print(round(n_active_cnes_at_saturation, digit=0))

print(round(n_active_cnes_at_saturation/NrCNEs, digit=4))

