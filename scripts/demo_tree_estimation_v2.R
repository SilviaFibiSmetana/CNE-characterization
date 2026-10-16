rm(list=ls())
args <- commandArgs ( TRUE )
if ( length ( args ) != 4 ){
	stop ( "Error. Check usage!" )
}
treefile <- args[1]
outputdir <- args[2]
orthtable <- args[3]
relspeciestable <- args[4]

if( ! file.exists(treefile) ){
  stop( paste(treefile, "does not exist!" ) )
}
if( ! file.exists(outputdir) ){
  stop( paste(outputdir, "does not exist!" ) )
}
if( ! file.exists(orthtable) ){
  stop( paste(orthtable, "does not exist!" ) )
}
if( ! file.exists(relspeciestable) ){
  stop( paste(relspeciestable, "does not exist!" ) )
}

library("dplyr")
library("tidyr")
library("ggplot2")
library("ape")
library("RRphylo")
library("data.table")
library("Matrix") 
library("pheatmap")
library("aplot")
library("ggtree")
library("phytools")
library("parallel")
library("tibble")

#prepare output
output <- paste(outputdir, basename(orthtable), "", sep="")

############################################################ 
### READ TREE 
############################################################ 
tree <- read.tree(treefile) 
TIPS <- rep(0, 470) 
names(TIPS) <- tree$tip.label 
############################################################ 
### DIVERGENCE TIMES 
############################################################ 
# Homininae = African apes = 8 Ma # Hominoidea = Apes = 18 Ma # Catarrhini = Old world monkeys + apes = 25 Ma # Primates = Primates = 65 Ma # Euarchontoglires = Rodents, primates, etc. = 75 Ma # Boreoeutheria = = 85 Ma # Eutheria = Placental mammals = 90 Ma # Theria = Therian mammals = 160 Ma # Mammalia = Mammals = 180 Ma 
NODES <- c( 180, 160, 90, 85, 75, 65, 25, 18, 8 ) 
names(NODES) <- c( 472, 473, 474, 475, 476, 479, 482, 483, 484 ) 
tree2 <- scaleTree( tree, tip.ages = TIPS, node.ages = NODES ) 
############################################################ 
### READ ORTHOLOG INFORMATION 
############################################################ 
cat("\n") 
cat("Reading ortholog table:\n") 
cat(orthtable) 
cat("\n\n") 
orth <- read.table( orthtable, sep = "\t", header = FALSE ) 
cat("Ortholog table dimensions:\n") 
print(dim(orth)) 
cat("\nFirst rows:\n") 
print(head(orth)) 
############################################################ 
### EXTRACT SPECIES AND CNE SETS 
############################################################ 
species <- readLines(relspeciestable)
species <- species[species != ""]
cat("Number of species in supplied species list:", length(species), "\n")
#species <- as.vector(unique(orth$V1)) 
sets <- as.vector(unique(orth$V2)) 
cat("\nNumber of species:") 
cat(length(species)) 
cat("\n") 
cat("Number of CNE sets:") 
cat(length(sets)) 
cat("\n") 
cat("\nCNE set IDs:\n") 
print(sets) 
############################################################ 
### GENERATE SPECIES x CNE MATRIX 
############################################################ 
zero_matrix <- matrix(
  0,
  nrow = length(species),
  ncol = length(sets)
)

df <- as.data.frame(zero_matrix)

colnames(df) <- as.character(sets)
rownames(df) <- species
#
############################################################ 
### CONVERT TO DATA.TABLE 
############################################################ 
orth_dt <- as.data.table(orth)

unique_V2 <- unique(orth_dt$V2)

dfsp <- sparseMatrix(
  i = integer(0),
  j = integer(0),
  dims = c(length(species), length(unique_V2))
)

rownames(dfsp) <- species

row_indices <- setNames(
  seq_along(species),
  species
)

for (i in seq_along(unique_V2)) {
  
  cne_id <- unique_V2[i]
  
  sub <- orth_dt[V2 == cne_id, V1]
  
  matched_species <- intersect(sub, species)
  
  if (length(matched_species) > 0) {
    
    matched_rows <- row_indices[matched_species]
    
    dfsp[matched_rows, i] <- 1
  }
}

df_test <- as.data.frame(as.matrix(dfsp))

colnames(df_test) <- as.character(unique_V2)
rownames(df_test) <- species
#orth_dt <- as.data.table(orth) 
#cat("\nFirst rows of ortholog data:\n") 
#print(head(orth_dt)) 
############################################################ 
### GET UNIQUE CNE SET IDs 
############################################################ 
#unique_V2 <- unique(orth_dt$V2) 
############################################################ 
### CREATE SPARSE MATRIX 
############################################################ 
#dfsp <- sparseMatrix( i = integer(0), j = integer(0), dims = c( nrow(df), length(unique_V2) ) ) 
#rownames(dfsp) <- rownames(df) 
############################################################ 
### SPECIES LOOKUP 
############################################################ 
#row_indices <- setNames( seq_len(nrow(dfsp)), rownames(dfsp) ) 
############################################################ 
### FILL SPECIES x CNE MATRIX 
############################################################ 
#cat("\n") 
#cat("Filling species x CNE matrix...\n") 
#for (i in seq_along(unique_V2)) { 
#  cne_id <- unique_V2[i] 
#  sub <- orth_dt[ V2 == cne_id, V1 ] 
#  matched_rows <- row_indices[ intersect( sub, rownames(dfsp) ) ] 
#  if (length(matched_rows) > 0) { 
#    dfsp[matched_rows, i] <- 1 
#    } 
#} 
#head(dfsp)
############################################################ 
### CONVERT TO REGULAR DATA FRAME 
############################################################ 
#df_test <- as.data.frame( as.matrix(dfsp) ) 
# Make sure columns have the actual CNE IDs 
#colnames(df_test) <- as.character(unique_V2) 
#rownames(df_test) <- rownames(dfsp) 
#head(df_test)
############################################################ 
### CHECK MATRIX 
############################################################ 
cat("\n") 
cat("Species x CNE matrix dimensions:\n") 
print(dim(df_test)) 
cat("\nNumber of CNE sets:") 
cat(ncol(df_test)) 
cat("\n") 
if (ncol(df_test) != length(unique(sets))) { 
  warning( "Number of columns in df_test differs from number of unique sets." ) } 
############################################################ 
### SPECIES x SPECIES MATRIX 
############################################################ 
cat("\n") 
cat("Generating species x species matrix...\n") 
specspec <- data.frame( matrix( nrow = length(species), ncol = length(species) ) ) 
rownames(specspec) <- species 
colnames(specspec) <- species 
specspec[,] <- as.matrix(df_test) %*% t(as.matrix(df_test)) 
cat("\nSpecies x species matrix dimensions:\n") 
print(dim(specspec)) 
############################################################ 
### REMOVE SPECIES NOT PRESENT IN ANALYSIS 
############################################################ 
all_tree_species <- tree2$tip.label

drop_tree_tips <- setdiff(all_tree_species, species)

tree3 <- drop.tip(tree2, drop_tree_tips)

length(tree3$tip.label)

setdiff(species, tree3$tip.label) #should give character(0)

cat("Species in supplied list:", length(species), "\n")
cat("Species in tree:", length(tree3$tip.label), "\n")

missing_from_tree <- setdiff(species, tree3$tip.label)

if (length(missing_from_tree) > 0) {
  stop(
    "Some species from species list are missing from the tree: ",
    paste(missing_from_tree, collapse=", ")
  )
}
#all_tree_species <- tree2$tip.label 
#drop_tree_tips <- setdiff( all_tree_species, species ) 
#tree3 <- drop.tip( tree2, drop_tree_tips ) 
############################################################ 
### TREE + HEATMAP OUTPUT 
############################################################ 
specspec1 <- specspec[ tree3$tip.label, tree3$tip.label ] 
write.table( specspec1, file = paste( output, "_speciesvsspecies.txt", sep = "" ), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE ) 
############################################################ 
### PREPARE TREE FOR ESTIMATION 
############################################################ 
set_tree <- tree3 
# Collapse singleton nodes 
set_tree <- collapse.singles(set_tree) 
cat("\n") 
cat("Number of species in tree:") 
cat(length(set_tree$tip.label)) 
cat("\n") 
############################################################ 
### INTERNAL NODES 
############################################################ 
internal_nodes <- unique( set_tree$edge[, 1] ) 
cat("\nNumber of internal nodes:") 
cat(length(internal_nodes)) 
cat("\n") 
############################################################ 
### PARALLEL ANCESTRAL ESTIMATION FUNCTION 
############################################################ 
process_set <- function( i, df_test, set_tree ) { 
  cat( "Processing CNE", 
       colnames(df_test)[i], "\n" ) 
  ########################################################## 
  ### GET CURRENT CNE 
  ########################################################## 
  set_binary <- df_test[, i] 
  ########################################################## 
  ### SPECIAL CASE: ### CNE PRESENT IN ALL 201 SPECIES 
  ########################################################## 
  if (sum(set_binary) == 201) { 
    res <- as.numeric( set_binary[-1] ) 
  ########################################################## 
    ### NORMAL CASE: ### ANCESTRAL CHARACTER ESTIMATION 
    ########################################################## 
    } else { 
      names(set_binary) <- rownames(df_test) 
      ######################################################## 
      ### MATCH SPECIES ORDER TO TREE 
      ######################################################## 
      set_binary <- set_binary[ set_tree$tip.label ] 
      ######################################################## 
      ### ACE 
      ######################################################## 
      reconstruction <- ace( set_binary, set_tree, type = "discrete" ) 
      ######################################################## 
      ### MOST LIKELY ANCESTRAL STATE 
      ######################################################## 
      res <- apply( reconstruction$lik.anc, 1, FUN = function(x) { which( x == max(x) ) - 1 } ) } 
  ########################################################## 
  ### RETURN 
  ########################################################## 
  return(res) } 
############################################################ 
### RUN PARALLEL ESTIMATION 
############################################################ 
num_cores <- detectCores() - 1 
cat("\n") 
cat("====================================================\n") 
cat("STARTING PARALLEL ANCESTRAL ESTIMATION\n") 
cat("====================================================\n") 
cat("CNE sets:", ncol(df_test), "\n") 
cat("Internal nodes:", length(internal_nodes), "\n") 
cat("Using cores:", num_cores, "\n") 
cat("====================================================\n\n") 
############################################################ 
### IMPORTANT: ### The input itself contains only 167 CNEs. 
### 
### Therefore: 
### - no batch calculation 
### - no 100000 offsets 
### - no artificial 1000-CNE chunks 
### 
### We simply parallelize over all CNEs. 
############################################################ 

process_set <- function(i, df_test, set_tree) {
  
  myset <- df_test[, i, drop = FALSE]
  set_binary <- myset[, 1]
  
  if (sum(set_binary) == 201) {
    return(as.numeric(set_binary[-1]))
  }
  
  names(set_binary) <- rownames(myset)
  set_binary <- set_binary[set_tree$tip.label]
  
  reconstruction <- ace(
    set_binary,
    set_tree,
    type = "discrete"
  )
  
  res <- apply(
    reconstruction$lik.anc,
    1,
    function(x) which.max(x) - 1
  )
  
  return(res)
}



results <- mclapply( seq_len(ncol(df_test)), 
                     process_set, 
                     df_test = df_test, 
                     set_tree = set_tree, 
                     mc.cores = num_cores ) 
############################################################ 
### COMBINE RESULTS 
############################################################ 
cat("\n") 
cat("Combining results...\n") 
resdfall <- as.data.frame( do.call( cbind, results ) ) 
############################################################ 
### ADD NAMES 
############################################################ 
rownames(resdfall) <- internal_nodes 
colnames(resdfall) <- colnames(df_test) 
############################################################ 
### CHECK FINAL RESULT 
############################################################ 
cat("\n") 
cat("====================================================\n") 
cat("FINAL RESULT\n") 
cat("====================================================\n") 
cat("Dimensions:\n") 
print(dim(resdfall)) 
cat("\nCNE sets processed:") 
cat(ncol(resdfall)) 
cat("\n") 
cat("Internal nodes:") 
cat(nrow(resdfall)) 
cat("\n") 
############################################################ 
### PREPARE OUTPUT 
############################################################ 
resdfallout <- cbind( RowNames = rownames(resdfall), resdfall ) 
############################################################ 
### WRITE TXT 
############################################################ 
outfile_txt <- paste( output, "_nodeestimations.txt", sep = "" ) 
cat("\nWriting output:\n") 
cat(outfile_txt) 
cat("\n") 
write.table( resdfallout, file = outfile_txt, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE ) 
############################################################ 
### SAVE RDATA 
############################################################ 
outfile_rdata <- paste( output, "_nodeestimations.RData", sep = "" ) 
cat("\nSaving RData:\n") 
cat(outfile_rdata) 
cat("\n") 
save( resdfallout, file = outfile_rdata ) 
############################################################ 
### DONE 
############################################################ 
cat("\n") 
cat("====================================================\n") 
cat("DONE\n") 
cat("====================================================\n") 
cat( "Successfully processed", ncol(resdfall), "CNE sets.\n" ) 
cat( "Estimated", nrow(resdfall), "internal nodes.\n" ) 
cat("====================================================\n")

