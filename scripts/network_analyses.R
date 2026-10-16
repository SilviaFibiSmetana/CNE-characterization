rm(list=ls())
args <- commandArgs(TRUE)
if (length (args) != 9){
	stop ("Error. Check usage!")
}
ACTIVITY_FILE <- args[1]
CLADES_FILE <- args[2]

PRIMATE_PROMOTER_FILE <- args[3]
MAMMAL_PROMOTER_FILE <- args[4]

INSOMNIA_FILE <- args[5]
SCHIZOPHRENIA_FILE <- args[6]
NEUROTICISM_FILE <- args[7]
NEUROTICISM_MEASUREMENT_FILE <- args[8]

OUTPUT_DIR <- args[9]

if( !file.exists(ACTIVITY_FILE) ){
  stop( paste(ACTIVITY_FILE, "does not exist!" ) )
}
if( !file.exists(CLADES_FILE) ){
  stop( paste(CLADES_FILE, "does not exist!" ) )
}

if( !file.exists(PRIMATE_PROMOTER_FILE) ){
  stop( paste(PRIMATE_PROMOTER_FILE, "does not exist!" ) )
}
if( !file.exists(MAMMAL_PROMOTER_FILE) ){
  stop( paste(MAMMAL_PROMOTER_FILE, "does not exist!" ) )
}

if( !file.exists(SCHIZOPHRENIA_FILE) ){
  stop( paste(SCHIZOPHRENIA_FILE, "does not exist!" ) )
}
if( !file.exists(NEUROTICISM_FILE) ){
  stop( paste(NEUROTICISM_FILE, "does not exist!" ) )
}
if( !file.exists(NEUROTICISM_MEASUREMENT_FILE) ){
  stop( paste(NEUROTICISM_MEASUREMENT_FILE, "does not exist!" ) )
}
if( !file.exists(OUTPUT_DIR) ){
  stop( paste(OUTPUT_DIR, "does not exist!" ) )
}

#DATA_DIR <- "./data"
#ACTIVITY_FILE <- file.path(DATA_DIR, "H3K27ac_DNase_CNE_activation.tsv")
#CLADES_FILE <- file.path(DATA_DIR, "summary_CNEs_in_clades_publication.txt")
#PRIMATE_PROMOTER_FILE <- file.path(DATA_DIR, "Primates_estimated_CNEs_small_prom_overlap")
#MAMMAL_PROMOTER_FILE <- file.path(DATA_DIR, "Mammals_estimated_CNEs_small_prom_overlap")
#INSOMNIA_FILE <- file.path(DATA_DIR, "all_estimated_CNEs_Insomnia_new_intersected_filtered")
#SCHIZOPHRENIA_FILE <- file.path(DATA_DIR, "all_estimated_CNEs_Schizophrenia_new_intersected_filtered")
#NEUROTICISM_FILE <- file.path(DATA_DIR, "all_estimated_CNEs_Neuroticism_new_intersected_filtered")
#NEUROTICISM_MEASUREMENT_FILE <- file.path(DATA_DIR, "all_estimated_CNEs_Neuroticismmeasurement_new_intersected_filtered")



# -----------------------------------------------------------------------------
# Paths and analysis parameters
# -----------------------------------------------------------------------------

output <- file.path(OUTPUT_DIR, "network_analysis.R")

set.seed(10)

N_PERCOLATION_LEVELS <- 100
N_PERCOLATION_NULL <- 100
N_SUBNETWORK_NULL <- 500
WEIGHT_THRESHOLD <- 0.38  

# -----------------------------------------------------------------------------
# Required libraries
# -----------------------------------------------------------------------------

library("Matrix")
library("dplyr")
library("tidyr")
library("tibble")
library("ggplot2")
library("ggthemes")
library("igraph")
library("vegan")
library("ggpubr") # P-values on the figure
library("rstatix") # tidy p-values
library("introdataviz") # split violins
library("patchwork")
library("tidygraph")
library("ggraph")
library("withr") # with_seed

# -----------------------------------------------------------------------------
# Visualization
# -----------------------------------------------------------------------------

my_clean_theme <- theme_clean() +
  theme(legend.position = "bottom",
    plot.background = element_blank(),
    panel.border = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    axis.title.x = element_blank()
  )

my_conservation_colors <- c("Mammals,Primates" = "#B6357F", "Primates" = "#7F7A43")

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

read_cne_coordinates <- function(path) {
  x <- read.delim(path, header = FALSE)
  unique(paste(paste(x[, 1], x[, 2], sep = ":"), x[, 3], sep = "-"))
}

ComputeNodeNetworkMetrics <- function(g){
	# 1. Degree
	deg <- degree(g)
	deg_stats <- c(mean = mean(deg), sd = sd(deg), median = median(deg))

	# 2. Strength
	str_val <- strength(g, weights = E(g)$weight)
	str_stats <- c(mean = mean(str_val), sd = sd(str_val), median = median(str_val))

	# 3. Normalized Weighted Betweenness
	bet <- betweenness(g, weights = E(g)$weight, normalized = TRUE)
	bet_stats <- c(mean = mean(bet), sd = sd(bet), median = median(bet))

	# 4. Weighted PageRank
	pr <- page_rank(g, weights = E(g)$weight)$vector
	pr_stats <- c(mean = mean(pr), sd = sd(pr), median = median(pr))

	# Combine results into a data frame for easy viewing
	metrics_summary <- data.frame(
	  Metric = c("Degree", "Strength", "Betweenness (Norm/W)", "PageRank (W)"),
	  Mean = c(deg_stats["mean"], str_stats["mean"], bet_stats["mean"], pr_stats["mean"]),
	  SD = c(deg_stats["sd"], str_stats["sd"], bet_stats["sd"], pr_stats["sd"]),
	  Median = c(deg_stats["median"], str_stats["median"], bet_stats["median"], pr_stats["median"])
	)

	return(metrics_summary)
}

ComputeGlobalNetworkMetrics <- function(g){
	# 1. Density
	dens <- edge_density(g)

	# 2. Transitivity (Global)
	trans <- transitivity(g, type = "global")

	# 3. Von Neumann Entropy (via Laplacian spectrum)
	# Obtain the unweighted or weighted unnormalized Laplacian matrix
	L <- laplacian_matrix(g, weights = E(g)$weight, normalization = "unnormalized", sparse = FALSE)
	# Compute eigenvalues
	evals <- eigen(L, symmetric = TRUE, only.values = TRUE)$values
	# Filter out near-zero eigenvalues to avoid log(0)
	evals <- evals[evals > 1e-10]
	# Normalize eigenvalues to form a probability distribution (density matrix trace = n or 1)
	# For standard graph von Neumann entropy, eigenvalues are scaled by sum of eigenvalues (or network volume)
	p <- evals / sum(evals)
	vn_entropy <- -sum(p * log(p))

	metrics_summary <- data.frame(Density = dens, Transitivity = trans, VonNeumannEntropy = vn_entropy)

	return(metrics_summary)
}

ComparePhenotypeCentrality <- function(filtered_g, network_cnes){
	filtered_g.df <- data.frame(CNE=names(V(filtered_g)))
	filtered_g.df <- filtered_g.df %>% left_join(network_cnes, by="CNE")
	#print(dim(filtered_g.df))

	phenotypes <- levels(factor(filtered_g.df$Phenotype))
	subnetworks <- list()
	induced_subnetworks_metrics <- list()
	node_metrics <- tibble()
	for (p in phenotypes){
		#print(p)
		nodes_to_keep <- filtered_g.df %>% filter(Phenotype==p) %>% select(CNE) %>% unique() %>% pull()
		sub_g <- induced_subgraph(filtered_g, nodes_to_keep)
		subnetworks[[p]] <- sub_g

		# Compute node-level metrics
		N <- vcount(sub_g)
		E <- ecount(sub_g)
		#print(N)
		#print(E)
		Degree <- degree(sub_g, mode = "all")/(N-1)
		Strength <- strength(sub_g, mode = "all", weights = E(sub_g)$weight)/sum(E(sub_g)$weight)
		PageRank_Centrality <- page_rank(sub_g, weights = E(sub_g)$weight)$vector
		Node_Betweenness <- betweenness(sub_g, normalized=TRUE)

		## Store individual node-level metrics in a table
		aux <- data.frame(cbind(Degree, Strength, PageRank_Centrality, Node_Betweenness))
		aux$CNE <- rownames(aux)
		rownames(aux) <- NULL
		aux$Phenotype <- p
		node_metrics <- bind_rows(node_metrics, aux)

		wtc <- cluster_louvain(sub_g, weights = E(sub_g)$weight)
		modularity_value <- modularity(wtc)

		# Calculate averages and store in a single-row data frame
		induced_subnetworks_metrics[[p]] <- data.frame(
		    nr_nodes = N,
		    nr_edges = E,
		    # Network connectivity
		    edge_density = edge_density(sub_g),
		    transitivity = transitivity(sub_g),
		    modularity = modularity_value,
		    # Node metrics
		    avg_Degree = mean(Degree),
		    avg_Strength = mean(Strength),
		    avg_PageRank_Centrality = mean(PageRank_Centrality),
		    avg_betweenness = mean(Node_Betweenness)
		)
	}

	node_metrics <- node_metrics %>% 
		left_join(network_cnes, by=c("CNE", "Phenotype"), relationship="many-to-many")


	# Long Format
	node_metrics.long <- node_metrics %>%
		pivot_longer( c(Degree, Strength, PageRank_Centrality, Node_Betweenness), values_to="value", names_to="metric" ) %>%
		mutate(metric=gsub("_", " ", metric))


	# Compute P-values
	stat_pvals <- node_metrics.long %>%
	  group_by(metric, Phenotype) %>%
	  wilcox_test(value ~ clade) %>%
	  ungroup() %>%
	  adjust_pvalue(method = "BH") %>%
	  add_significance() %>%
	  add_y_position()

	p <- ggplot(node_metrics.long %>% filter(isPromoter=="Enhancer"), 
	       aes(x=Phenotype, y=value, fill=clade)) +
	       geom_split_violin(alpha = 1, trim = FALSE, aes(col=clade, group=interaction(clade, Phenotype))) +
	     geom_boxplot(
	    aes(fill=clade),
	    width = 0.1,
	    alpha = 0.5,
	    fatten = NULL,
	    show.legend = FALSE,
	    position = position_dodge(0.175) # Matches your stat_summary dodge
	  ) +
	  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
		       position = position_dodge(.175)) +
	  labs(x="Degree Centrality") +
	#  scale_x_discrete(name = "Condition", labels = c("Non-word", "Word")) +
	#  scale_y_continuous(name = "Reaction time (ms)",
	#                     breaks = seq(200, 800, 100), 
	#                     limits = c(200, 800) ) +
	 # scale_fill_brewer(palette = "Dark2", name = "Language group") +
	  scale_fill_manual(values = c("#B6357F", "#7F7A43")) +
	  scale_color_manual(values = c("#B6357F", "#7F7A43")) +
	  # Add the P-values
	#  stat_pvalue_manual(stat_pvals, 
	#		     label = "p.adj.signif", 
	#		     tip.length = 0.01) +
	  geom_text(
	  data = stat_pvals,
	  aes(
	    x = Phenotype,
	    y = y.position,
	    label = p.adj.signif    # The significance stars (*, **, ns)
	  ),
	  inherit.aes = FALSE,      # Keeps the text from trying to split like the violins
	#  vjust = 0.3,             # Nudges the stars just above the violin tip
	  size = 4 
	) +
	  my_clean_theme +
	  theme(legend.title = element_blank()) +
	  facet_wrap(~metric, scales = "free_y", nrow=1)
	output_pdf_file <- paste(output, "_CentralityVsConservation.pdf", sep="")
	pdf(output_pdf_file, width=12, height=3)
	print(p)
	dev.off()
}

network_prunning <- function(g, output){
	# 0. Plot Histogram of Edge Weights
	edge_df <- data.frame(Weight = E(g)$weight)

	p1 <- ggplot(edge_df, aes(x = Weight)) +
	  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white", alpha = 0.8) +
	  theme_clean() +
	  labs(
	#    title = "Distribution of Network Edge Weights",
	    x = "Edge Weight",
	    y = "Count"
	  )

	# 1. Define your threshold range (e.g., from 0 to 1 by steps of 0.05)
	# Adjust these numbers based on your actual edge weight distribution
	thresholds <- seq(0, 1, by = 0.01)

	# 2. Loop through thresholds and calculate network metrics
	metrics_history <- lapply(thresholds, function(t) {
	  # Prune edges below the threshold
	  filtered_g <- delete_edges(g, E(g)[weight < t])
	  
	  # Calculate Modularity (using Louvain community detection)
	  # Handle empty graphs or isolated nodes gracefully
	  mod_value <- if(vcount(filtered_g) > 0 && ecount(filtered_g) > 0) {
	    modularity(cluster_louvain(filtered_g))
	  } else {
	    0
	  }
	  
	  # Return a single row data frame for this threshold
	  data.frame(
	    Threshold    = t,
	    Density      = edge_density(filtered_g),
	    Transitivity = transitivity(filtered_g, type = "global"),
	    Modularity   = mod_value
	  )
	}) %>% bind_rows()

	# 3. Reshape the data for ggplot2 (long format)
	metrics_long <- metrics_history %>%
	  pivot_longer(cols = c(Density, Transitivity, Modularity),
		       names_to = "Metric",
		       values_to = "Value")

	# 4. Plot the evolution of connectivity metrics
	p2 <- ggplot(metrics_long, aes(x = Threshold, y = Value)) +
	  geom_line(linewidth = 1, color="red") +
	#  geom_point(size = 2, color="red") +
	  facet_wrap(~Metric, scales = "free_y", nrow = 1) +
	  theme_clean() +
	  labs(
	#    title = "Network Connectivity Evolution During Pruning",
	    x = "Edge Weight Pruning Threshold",
	    y = "Metric Value"
	  ) +
	  theme(
	    legend.position = "none",
	    strip.text = element_text(size = 12)
	  )

	# Display the combined layout
	p <- (p1 | p2) + 
	  plot_layout(widths = c(0.5, 1.5)) &
	  theme(
	    panel.border     = element_blank(),  # Removes plot panel frame
	    plot.background  = element_blank(),  # Removes outer background frame
	    axis.line        = element_line(color = "black") # Optional: adds clean axis lines instead
	  )
	  #plot_annotation(title = "Network Threshold Optimization Analysis")

	output_pdf_file <- paste(output, "_NetworkPrunning.pdf", sep="")
	pdf(output_pdf_file, width=10, height=3)
	print(p)
	dev.off()
}

plot_subnetwork <- function(my_network, coords_df, title, output) {

  my_conservation_colors <- c("Mammals,Primates" = "#B6357F", "Primates" = "#7F7A43")

  # Convert subnetwork to tidygraph
  sub_tg <- as_tbl_graph(my_network)

  # Filter the master coordinates to only include nodes present in this subnetwork
  # and sort them to match the exact order of nodes in sub_tg
  node_names <- with_graph(sub_tg, .N()$name)
  matched_coords <- coords_df[match(node_names, coords_df$name), c("x", "y")]

  # Pass the matched matrix/data.frame directly into the layout argument
  p <- ggraph(sub_tg, layout = matched_coords) +
    geom_edge_link(aes(width = weight), alpha = 0.2, color = "grey40") +
    scale_edge_width(range = c(0.2, 2)) +
    geom_node_point(aes(color = conservation, fill = conservation, shape = trait), size=3) +

    # Outline for promoters
    geom_node_point(
    data = function(x) x[!is.na(x$GeneSymbol), ], # Filters data to ONLY valid genes
    shape = 3,
    size = 1,
    color = "black",                               # High-contrast accent color for the ring
    fill = "transparent", # Hollow interior so the core node is visible
    stroke = 0.5 # Thick outline border ring
  ) +
    geom_node_text(aes(label = GeneSymbol), repel = TRUE, family = "sans", fontface = "italic", size = 3) +

    scale_size_continuous(range = c(2, 8)) +
    scale_fill_manual(values = my_conservation_colors) +
    scale_color_manual(values = my_conservation_colors) +
    scale_shape_manual(values = c(21, 8, 9, 25, 7, 22)) +
    theme_graph(base_family = 'sans')

    output_pdf_file <- paste(output, "_NetworkGraph.pdf", sep="")
    pdf(output_pdf_file, width=10, height=10)
    print(p)
    dev.off()
}






von_neumann_entropy <- function(g, weight_attr = "weight", tol = 1e-12) {
  stopifnot(!is_directed(g))
  if (vcount(g) == 0 || ecount(g) == 0) return(0)

  w <- edge_attr(g, weight_attr)
  if (is.null(w)) stop(sprintf("Edge attribute '%s' not found.", weight_attr))
  if (any(w < 0, na.rm = TRUE)) stop("Von Neumann entropy requires non-negative weights.")

  keep <- is.finite(w) & w > 0
  if (!all(keep)) g <- delete_edges(g, E(g)[!keep])
  if (ecount(g) == 0) return(0)

  w <- edge_attr(g, weight_attr)
  L <- laplacian_matrix(g, weights = w, normalized = FALSE, sparse = FALSE)
  trL <- sum(diag(L))
  if (trL <= 0) return(0)

  mu <- eigen(L, symmetric = TRUE, only.values = TRUE)$values
  mu[mu < 0 & mu > -sqrt(tol)] <- 0
  mu <- mu[mu > tol]
  if (length(mu) == 0) return(0)

  lambda <- mu / trL
  -sum(lambda * log(lambda))
}

node_metrics <- function(g) {
  w <- E(g)$weight
  if (length(w) > 0 && any(w <= 0 | !is.finite(w))) {
    stop("All retained network weights must be finite and > 0.")
  }

  # Convert similarity weights to distances for weighted betweenness.
  distances <- if (length(w) > 0) 1 / w else numeric(0)

  data.frame(
    node = V(g)$name,
    node_type = V(g)$node_type,
    conservation = V(g)$conservation,
    degree = degree(g, mode = "all"),
    strength = strength(g, mode = "all", weights = w),
    pagerank = page_rank(g, directed = FALSE, weights = w)$vector,
    betweenness = betweenness(
      g,
      directed = FALSE,
      weights = distances,
      normalized = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

network_metrics <- function(g) {
  w <- E(g)$weight
  comp <- components(g)
  louvain <- if (ecount(g) > 0 && vcount(g) > 1) {
    cluster_louvain(g, weights = w)
  } else {
    NULL
  }

  data.frame(
    nodes = vcount(g),
    edges = ecount(g),
    density = edge_density(g, loops = FALSE),
    mean_degree = mean(degree(g, mode = "all")),
    transitivity = if (ecount(g) > 0) transitivity(g, type = "global") else NA_real_,
    components = comp$no,
    largest_component_share = if (vcount(g) > 0) max(comp$csize) / vcount(g) else NA_real_,
    modularity = if (!is.null(louvain)) modularity(louvain) else NA_real_,
    von_neumann_entropy = von_neumann_entropy(g),
    stringsAsFactors = FALSE
  )
}

percolation_metrics <- function(g) {
  n <- vcount(g)
  deg <- degree(g)
  comp <- components(g)
  lcc_vertices <- which(comp$membership == which.max(comp$csize))
  lcc <- induced_subgraph(g, vids = lcc_vertices)
  q <- if (ecount(lcc) > 0 && vcount(lcc) > 1) {
    modularity(cluster_louvain(lcc, weights = E(lcc)$weight))
  } else {
    NA_real_
  }

  data.frame(
    components = comp$no,
    isolates = sum(deg == 0),
    largest_component_share = if (n > 0) max(comp$csize) / n else NA_real_,
    modularity = q
  )
}

weight_prune_curve <- function(g, steps = 100) {
  m <- ecount(g)
  if (m == 0) stop("Percolation analysis requires at least one edge.")

  edge_order <- order(E(g)$weight, decreasing = FALSE)
  ordered_weights <- E(g)$weight[edge_order]
  fractions <- seq(0, 1, length.out = steps)

  bind_rows(lapply(fractions, function(f) {
    k <- floor(f * m)
    gk <- if (k == 0) g else delete_edges(g, E(g)[edge_order[seq_len(k)]])
    metrics <- percolation_metrics(gk)
    metrics$frac_removed <- f
    metrics$edges_removed <- k
    metrics$weight_threshold_removed <- if (k == 0) NA_real_ else ordered_weights[k]
    metrics
  }))
}

random_prune_curve <- function(g, steps = 100, n_sim = 100, seed = 10) {
  m <- ecount(g)
  if (m == 0) stop("Percolation analysis requires at least one edge.")

  fractions <- seq(0, 1, length.out = steps)
  k_removed <- floor(fractions * m)
  set.seed(seed)

  all_null <- bind_rows(lapply(seq_len(n_sim), function(sim) {
    edge_order <- sample.int(m)
    bind_rows(lapply(seq_along(k_removed), function(i) {
      k <- k_removed[i]
      gk <- if (k == 0) g else delete_edges(g, E(g)[edge_order[seq_len(k)]])
      out <- percolation_metrics(gk)
      out$sim <- sim
      out$frac_removed <- fractions[i]
      out$edges_removed <- k
      out
    }))
  }))

  all_null %>%
    group_by(frac_removed, edges_removed) %>%
    summarise(
      components_mean = mean(components, na.rm = TRUE),
      components_q05 = quantile(components, 0.05, na.rm = TRUE),
      components_q95 = quantile(components, 0.95, na.rm = TRUE),
      isolates_mean = mean(isolates, na.rm = TRUE),
      isolates_q05 = quantile(isolates, 0.05, na.rm = TRUE),
      isolates_q95 = quantile(isolates, 0.95, na.rm = TRUE),
      largest_component_share_mean = mean(largest_component_share, na.rm = TRUE),
      largest_component_share_q05 = quantile(largest_component_share, 0.05, na.rm = TRUE),
      largest_component_share_q95 = quantile(largest_component_share, 0.95, na.rm = TRUE),
      modularity_mean = mean(modularity, na.rm = TRUE),
      modularity_q05 = quantile(modularity, 0.05, na.rm = TRUE),
      modularity_q95 = quantile(modularity, 0.95, na.rm = TRUE),
      .groups = "drop"
    )
}

sample_matched_nodes <- function(g, target_nodes) {
  target_type <- V(g)$node_type[match(target_nodes, V(g)$name)]
  need <- table(target_type)

  sampled <- unlist(lapply(names(need), function(type) {
    pool <- V(g)$name[V(g)$node_type == type]
    sample(pool, size = as.integer(need[[type]]), replace = FALSE)
  }), use.names = FALSE)

  sampled
}

subnetwork_descriptors <- function(g) {
  nm <- node_metrics(g)
  gm <- network_metrics(g)

  c(
    density = gm$density,
    transitivity = gm$transitivity,
    modularity = gm$modularity,
    von_neumann_entropy = gm$von_neumann_entropy,
    degree_mean = mean(nm$degree, na.rm = TRUE),
    strength_mean = mean(nm$strength, na.rm = TRUE),
    pagerank_mean = mean(nm$pagerank, na.rm = TRUE),
    betweenness_mean = mean(nm$betweenness, na.rm = TRUE)
  )
}

# Helper functions for phenotype subnetwork metrics and matched null sampling.
safe_val <- function(expr, default = NA_real_) {
  tryCatch(expr, error = function(e) default)
}

integrated_node_metrics <- function(g) {
  w <- igraph::edge_attr(g, "weight")
  if (is.null(w)) stop("Edge attribute 'weight' not found.")
  if (anyNA(w)) stop("Edge weights contain NA.")

  data.frame(
    node = V(g)$name,
    type = factor(V(g)$indicator),
    degree = igraph::degree(g, mode = "all"),
    strength = igraph::strength(g, mode = "all", weights = w),
    eigenvector_w = safe_val(
      igraph::eigen_centrality(g, weights = w)$vector
    ),
    pagerank_w = safe_val(
      igraph::page_rank(g, weights = w)$vector
    ),
    stringsAsFactors = FALSE
  )
}

sample_nodes_stratified <- function(
  g0,
  target_nodes,
  name_attr = "name",
  indicator_attr = "indicator",
  exact_indicator_levels = FALSE,
  group_map = NULL
) {
  vnames <- igraph::vertex_attr(g0, name_attr)
  if (is.null(vnames)) stop(sprintf("Vertex attribute '%s' not found.", name_attr))
  if (anyDuplicated(vnames)) stop(sprintf("Vertex attribute '%s' has duplicates.", name_attr))

  ind0 <- igraph::vertex_attr(g0, indicator_attr)
  if (is.null(ind0)) stop(sprintf("Vertex attribute '%s' not found.", indicator_attr))
  names(ind0) <- vnames

  if (any(!(target_nodes %in% vnames))) {
    missing <- setdiff(target_nodes, vnames)
    stop(sprintf(
      "Some target_nodes not found in g0 (first up to 10): %s",
      paste(head(missing, 10), collapse = ", ")
    ))
  }

  ind_t <- ind0[target_nodes]

  if (exact_indicator_levels) {
    grp0 <- as.character(ind0)
    grpt <- as.character(ind_t)
  } else {
    if (is.null(group_map)) {
      group_map <- function(x) {
        ifelse(
          x %in% c(1, 2), "enhancer",
          ifelse(x %in% c(3, 4), "promoter", NA_character_)
        )
      }
    }
    grp0 <- group_map(ind0)
    names(grp0) <- names(ind0)
    grpt <- group_map(ind_t)
    names(grpt) <- names(ind_t)
  }

  if (any(is.na(grpt))) stop("Some target nodes have indicator values that do not map to a group.")
  if (any(is.na(grp0))) stop("Some nodes in g0 have indicator values that do not map to a group.")

  need <- table(grpt)
  picked <- character(0)

  for (lab in names(need)) {
    pool <- names(grp0)[grp0 == lab]
    if (length(pool) < need[[lab]]) {
      stop(sprintf(
        "Not enough nodes in group '%s' to sample: need=%d, have=%d",
        lab, need[[lab]], length(pool)
      ))
    }
    picked <- c(
      picked,
      sample(pool, size = need[[lab]], replace = FALSE)
    )
  }

  picked
}

graph_modularity_louvain <- function(
  g,
  weight_attr = "weight",
  nstart = 20
) {
  stopifnot(!igraph::is_directed(g))

  w <- igraph::edge_attr(g, weight_attr)
  if (!is.null(w) && any(w <= 0, na.rm = TRUE)) {
    g <- igraph::delete_edges(g, E(g)[w <= 0])
  }

  if (igraph::ecount(g) == 0) {
    return(list(
      Q = NA_real_,
      membership = rep(NA_integer_, igraph::vcount(g))
    ))
  }

  best_Q <- -Inf
  best_mem <- NULL

  for (i in seq_len(nstart)) {
    cl <- igraph::cluster_louvain(
      g,
      weights = igraph::edge_attr(g, weight_attr)
    )
    Q <- igraph::modularity(cl)

    if (is.finite(Q) && Q > best_Q) {
      best_Q <- Q
      best_mem <- igraph::membership(cl)
    }
  }

  list(Q = best_Q, membership = best_mem)
}

integrated_subnetwork_significance <- function(
  g0,
  sub_nodes,
  n_sim = 500,
  weight_attr = "weight",
  name_attr = "name",
  indicator_attr = "indicator",
  exact_indicator_levels = FALSE,
#  metrics = c("degree", "strength", "eigenvector_w", "pagerank_w"),
  metrics = c("degree", "strength", "pagerank_w"),
  summaries = c("mean"),
  include_entropy = TRUE,
  entropy_normalized_laplacian = FALSE,
  include_modularity = TRUE,
  modularity_nstart = 20,
  alternative = c("greater", "less", "two.sided"),
  tol = 1e-12,
  p_adjust_method = c("none", "BH", "bonferroni", "holm")
) {
  alternative <- match.arg(alternative)
  p_adjust_method <- match.arg(p_adjust_method)

  vnames <- igraph::vertex_attr(g0, name_attr)
  if (is.null(vnames)) stop(sprintf("Vertex attribute '%s' not found in g0.", name_attr))

  if (any(!(sub_nodes %in% vnames))) {
    missing <- setdiff(sub_nodes, vnames)
    stop(sprintf(
      "Some sub_nodes not found in g0 (first up to 10): %s",
      paste(head(missing, 10), collapse = ", ")
    ))
  }

  g1 <- igraph::induced_subgraph(
    g0,
    vids = V(g0)[V(g0)$name %in% sub_nodes]
  )

  w1 <- igraph::edge_attr(g1, weight_attr)
  if (!is.null(w1) && any(w1 <= 0, na.rm = TRUE)) {
    g1 <- igraph::delete_edges(g1, E(g1)[w1 <= 0])
  }

  obs_node <- integrated_node_metrics(g1)

  obs <- list(
    n = igraph::vcount(g1),
    m = igraph::ecount(g1),
    edge_density = igraph::edge_density(g1, loops = FALSE)
  )

  for (m in metrics) {
    x <- obs_node[[m]]
    if ("mean" %in% summaries) {
      obs[[paste0(m, "_mean")]] <- mean(x, na.rm = TRUE)
    }
    if ("median" %in% summaries) {
      obs[[paste0(m, "_median")]] <- median(x, na.rm = TRUE)
    }
  }

  if (include_entropy) {
    obs$von_neumann_entropy <- von_neumann_entropy(
      g1,
      weight_attr = weight_attr,
      tol = tol
    )
  }

  if (include_modularity) {
    obs$modularity <- graph_modularity_louvain(
      g1,
      weight_attr = weight_attr,
      nstart = modularity_nstart
    )$Q
  }

  desc_names <- setdiff(names(obs), c("n", "m"))
  null_mat <- matrix(
    NA_real_,
    nrow = n_sim,
    ncol = length(desc_names),
    dimnames = list(NULL, desc_names)
  )

  for (s in seq_len(n_sim)) {
    rnd_nodes <- sample_nodes_stratified(
      g0,
      target_nodes = sub_nodes,
      name_attr = name_attr,
      indicator_attr = indicator_attr,
      exact_indicator_levels = exact_indicator_levels
    )

    gs <- igraph::induced_subgraph(
      g0,
      vids = V(g0)[V(g0)$name %in% rnd_nodes]
    )

    ws <- igraph::edge_attr(gs, weight_attr)
    if (!is.null(ws) && any(ws <= 0, na.rm = TRUE)) {
      gs <- igraph::delete_edges(gs, E(gs)[ws <= 0])
    }

    null_mat[s, "edge_density"] <- igraph::edge_density(gs, loops = FALSE)

    null_node <- integrated_node_metrics(gs)

    for (m in metrics) {
      x <- null_node[[m]]
      if ("mean" %in% summaries) {
        null_mat[s, paste0(m, "_mean")] <- mean(x, na.rm = TRUE)
      }
      if ("median" %in% summaries) {
        null_mat[s, paste0(m, "_median")] <- median(x, na.rm = TRUE)
      }
    }

    if (include_entropy) {
      null_mat[s, "von_neumann_entropy"] <- von_neumann_entropy(
        gs,
        weight_attr = weight_attr,
        tol = tol
      )
    }

    if (include_modularity) {
      null_mat[s, "modularity"] <- graph_modularity_louvain(
        gs,
        weight_attr = weight_attr,
        nstart = modularity_nstart
      )$Q
    }
  }

  emp_p <- function(sim_vec, obs_val, alternative) {
    sim_vec <- sim_vec[is.finite(sim_vec)]
    if (length(sim_vec) == 0 || !is.finite(obs_val)) return(NA_real_)

    if (alternative == "greater") {
      (sum(sim_vec >= obs_val) + 1) / (length(sim_vec) + 1)
    } else if (alternative == "less") {
      (sum(sim_vec <= obs_val) + 1) / (length(sim_vec) + 1)
    } else {
      mu <- mean(sim_vec)
      dev_obs <- abs(obs_val - mu)
      (sum(abs(sim_vec - mu) >= dev_obs) + 1) / (length(sim_vec) + 1)
    }
  }

  p <- vapply(
    colnames(null_mat),
    function(nm) emp_p(null_mat[, nm], obs[[nm]], alternative),
    numeric(1)
  )

  p_adj <- if (p_adjust_method == "none") {
    rep(NA_real_, length(p))
  } else {
    p.adjust(p, method = p_adjust_method)
  }

  null_mean <- colMeans(null_mat, na.rm = TRUE)
  null_median <- apply(null_mat, 2, median, na.rm = TRUE)
  null_sd <- apply(null_mat, 2, sd, na.rm = TRUE)

  results <- data.frame(
    descriptor = colnames(null_mat),
    observed = unname(vapply(
      colnames(null_mat),
      function(nm) obs[[nm]],
      numeric(1)
    )),
    null_mean = unname(null_mean),
    null_median = unname(null_median),
    null_sd = unname(null_sd),
    p_value = unname(p),
    p_adj = unname(p_adj),
    stringsAsFactors = FALSE
  )

  results$z <- (
    results$observed - results$null_mean
  ) / results$null_sd

  list(
    g1 = g1,
    observed = obs,
    null = null_mat,
    results = results
  )
}

make_phenotype_graph <- function(g_exp, phenotype_nodes) {
  ends_mat <- igraph::ends(g_exp, E(g_exp))

  keep_edge <- (
    ends_mat[, 1] %in% phenotype_nodes &
    ends_mat[, 2] %in% phenotype_nodes
  )

  gp <- igraph::delete_edges(
    g_exp,
    E(g_exp)[!keep_edge]
  )

  gp <- igraph::delete_vertices(
    gp,
    V(gp)[igraph::degree(gp) == 0]
  )

  ord <- order(V(gp)$name)
  igraph::permute(gp, ord)
}


# -----------------------------------------------------------------------------
# 1. CNE activity and phenotype-associated node sets
# -----------------------------------------------------------------------------

# Retain CNEs active in at least one tissue.
activity <- read.delim(ACTIVITY_FILE, header = TRUE, row.names = 1, check.names = FALSE)
activity <- as.matrix(activity) * 1
storage.mode(activity) <- "integer"
activity <- Matrix(activity, sparse = TRUE)
activity <- activity[rowSums(activity != 0) > 0, , drop = FALSE]

insomnia <- read_cne_coordinates(INSOMNIA_FILE)
schizophrenia <- read_cne_coordinates(SCHIZOPHRENIA_FILE)
neuroticism <- unique(c(
  read_cne_coordinates(NEUROTICISM_FILE),
  read_cne_coordinates(NEUROTICISM_MEASUREMENT_FILE)
))


phenotype_nodes <- list(
  Insomnia = sort(intersect(insomnia, rownames(activity))),
  Schizophrenia = sort(intersect(schizophrenia, rownames(activity))),
  Neuroticism = sort(intersect(neuroticism, rownames(activity)))
)


all_disease_nodes <- unique(unlist(phenotype_nodes, use.names = FALSE))
activity_disease <- activity[all_disease_nodes, , drop = FALSE]

# Promoter/enhancer annotation retained because it is required for matched null
# sampling and for the enhancer-only mammal/primate comparison.
pprom <- read.delim(PRIMATE_PROMOTER_FILE, header = FALSE)
mprom <- read.delim(MAMMAL_PROMOTER_FILE, header = FALSE)
promoter_table <- bind_rows(pprom, mprom) %>%
	unite(aux, V1, V2, sep=":") %>% 
	unite(CNE, aux, V3, sep="-") %>%
	rename(GeneSymbol=V4) %>%
	filter(GeneSymbol!=".") %>% 
	unique()
print(dim(promoter_table))

promoter_nodes <- unique(promoter_table$CNE)

clades <- read.delim(CLADES_FILE, header = TRUE, row.names = 1, check.names = FALSE)
clades[is.na(clades)] <- 0
mammal_nodes <- rownames(clades)[clades$Mammals == 1]
primate_nodes <- rownames(clades)[clades$Primates == 1]
primate_specific_nodes <- setdiff(primate_nodes, mammal_nodes)

# -----------------------------------------------------------------------------
# 2. Network CNE function and conservation
# -----------------------------------------------------------------------------

phenotype_nodes.long <- enframe(phenotype_nodes) %>% 
	unnest_longer(value) %>%
	rename(Phenotype=name, CNE=value)

clades.long <- clades %>% 
	rownames_to_column(var="CNE") %>% 
	pivot_longer(!CNE, names_to="clade", values_to="value") %>%
	filter(value!=0) %>%
	select(-value) %>%
	filter(clade %in% c("Mammals", "Primates"))
print(dim(clades.long))

network_cnes <- data.frame(CNE=rownames(activity_disease))
network_cnes <- network_cnes %>% 
	mutate(isPromoter=ifelse(CNE %in% promoter_nodes, "Promoter", "Enhancer")) %>%
	left_join(phenotype_nodes.long, by="CNE") %>%
	left_join(clades.long, by="CNE", relationship="many-to-many") %>%
	group_by(CNE, isPromoter, Phenotype) %>%
	summarise(clade=paste(sort(clade), collapse=",")) %>%
	left_join(promoter_table, by="CNE")
print(head(network_cnes))
print(dim(network_cnes))

print(network_cnes %>% group_by(Phenotype, isPromoter) %>% summarise(n_distinct(CNE)))

aux <- network_cnes %>%
        group_by(Phenotype, isPromoter, GeneSymbol, clade) %>% 
	summarise(NrCNEs=n_distinct(CNE)) %>% 
	mutate(FractionCNEs=NrCNEs/sum(NrCNEs)) %>%
        mutate(across(where(is.numeric), round, 3))
#print(aux)

p <- ggplot(aux,
       aes(x=isPromoter, y=NrCNEs, fill=clade)) +
        geom_col() +
        labs(y="Number of CNEs") +
        my_clean_theme +
        #scale_fill_npg() +
    scale_fill_manual(values = my_conservation_colors) +
    facet_wrap(~Phenotype) +
        theme(legend.title=element_blank())

output_pdf_file <- paste(output, "_NetworkCNEConservationByTrait.pdf", sep="")
pdf(output_pdf_file, width=9, height=5)
print(p)
dev.off()

# -----------------------------------------------------------------------------
# 3. CNE-based disease regulatory network construction
# -----------------------------------------------------------------------------

jaccard_dissimilarity <- vegan::vegdist(
  as.matrix(activity_disease),
  method = "jaccard",
  binary = TRUE
)

jaccard_similarity <- 1 - as.matrix(jaccard_dissimilarity)
diag(jaccard_similarity) <- 0
jaccard_similarity[!is.finite(jaccard_similarity)] <- 0

g_full <- graph_from_adjacency_matrix(
  jaccard_similarity,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

V(g_full)$node_type <- ifelse(V(g_full)$name %in% promoter_nodes, "promoter", "enhancer")
V(g_full)$conservation <- ifelse(
  V(g_full)$name %in% mammal_nodes,
  "mammal",
  ifelse(V(g_full)$name %in% primate_specific_nodes, "primate_specific", "other")
)

# Node classes: 1/2 = primate-specific/mammalian enhancers
# 3/4 = primate-specific/mammalian promoters.
V(g_full)$indicator <- (
  as.integer(V(g_full)$name %in% promoter_nodes) * 2L
) + as.integer(V(g_full)$name %in% mammal_nodes) + 1L


write("# Number of nodes and edges in the full network:", file="")
print(vcount(g_full))
print(ecount(g_full)) 

g_full_node_metrics <- ComputeNodeNetworkMetrics(g_full)

print(g_full_node_metrics %>% mutate(across(where(is.numeric), ~ round(., 4))))

g_full_global_metrics <- ComputeGlobalNetworkMetrics(g_full)
print(g_full_global_metrics %>% mutate(across(where(is.numeric), ~ round(., 4))))


# -----------------------------------------------------------------------------
# 4. Network robustness analysis and filtering
# -----------------------------------------------------------------------------

percolation_weighted <- weight_prune_curve(g_full, steps = N_PERCOLATION_LEVELS)
percolation_random <- random_prune_curve(
  g_full,
  steps = N_PERCOLATION_LEVELS,
  n_sim = N_PERCOLATION_NULL,
  seed = 10
)

percolation_long <- percolation_weighted %>%
  select(frac_removed, components, isolates, largest_component_share, modularity) %>%
  pivot_longer(-frac_removed, names_to = "metric", values_to = "observed") %>%
  left_join(
    percolation_random %>%
      select(
        frac_removed,
        components_mean, components_q05, components_q95,
        isolates_mean, isolates_q05, isolates_q95,
        largest_component_share_mean, largest_component_share_q05, largest_component_share_q95,
        modularity_mean, modularity_q05, modularity_q95
      ) %>%
      pivot_longer(
        -frac_removed,
        names_to = c("metric", ".value"),
        names_pattern = "(components|isolates|largest_component_share|modularity)_(mean|q05|q95)"
      ),
    by = c("frac_removed", "metric")
  )

p_percolation <- ggplot(percolation_long, aes(x = frac_removed)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.2) +
  geom_line(aes(y = mean), linetype = "dashed") +
  geom_line(aes(y = observed)) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(
    x = "Proportion of edges removed",
    y = NULL,
    title = "Network percolation: ordered low-weight removal vs random removal"
  ) +
  theme_bw()

ggsave(
  file.path(OUTPUT_DIR, "percolation_analysis.pdf"),
  plot = p_percolation,
  width = 9,
  height = 7
)


# Supplementary Figure
network_prunning(g_full, output)



# Apply the fixed similarity cutoff, retaining vertices that become isolated.
g <- igraph::delete_edges(
  g_full,
  igraph::E(g_full)[igraph::E(g_full)$weight < WEIGHT_THRESHOLD]
)

if (igraph::vcount(g) != igraph::vcount(g_full)) {
  stop("Thresholding changed the vertex count. Isolated vertices must be retained.")
}

message(sprintf(
  "# The pruned network has %d nodes and %d edges (Jaccard >= %.2f).",
  igraph::vcount(g), igraph::ecount(g), WEIGHT_THRESHOLD
))

# -----------------------------------------------------------------------------
# 5. Filtered network description
# -----------------------------------------------------------------------------

g_metrics <- ComputeNodeNetworkMetrics(g)

print(g_metrics %>% mutate(across(where(is.numeric), ~ round(., 4))))

# Table S6
write.csv(
  g_metrics %>% mutate(across(where(is.numeric), ~ round(., 4))),
  file.path(OUTPUT_DIR, "pruned_network_metrics_TableS6.csv"),
  row.names = FALSE
)

g_global_metrics <- ComputeGlobalNetworkMetrics(g)
print(g_global_metrics %>% mutate(across(where(is.numeric), ~ round(., 4))))

## Visualize Graph

# Prepare data for visualization
filtered_g <- g
## Add gene symbols
aux <- V(filtered_g)$name
node_annotation <- network_cnes %>% 
    group_by(CNE) %>% 
    mutate(Phenotype=paste(sort(Phenotype), collapse=",")) %>%
    ungroup() %>%
    select(CNE, GeneSymbol, clade, Phenotype) %>%
    unique() %>%
    filter(CNE %in% aux) %>%
    arrange(match(CNE, aux)); 
if( !all(aux == aux) ){
    stop("Check your data!")
}
# Gene Symbol
aux <- node_annotation$GeneSymbol
V(filtered_g)$GeneSymbol <- aux
# Conservation
aux <- node_annotation$clade
V(filtered_g)$conservation <- aux
# Phenotype
aux <- node_annotation$Phenotype
V(filtered_g)$trait <- aux
filtered_g.df <- as_tbl_graph(filtered_g)
#print(filtered_g.df)

# Generate a master coordinate data frame from the parent layout
#parent_layout <- create_layout(as_tbl_graph(filtered_g), layout="fr")
with_seed(42, {
	parent_layout <- create_layout(as_tbl_graph(filtered_g), layout="stress")
	master_coords <- data.frame(
	  name = parent_layout$name,
	  x = parent_layout$x,
	  y = parent_layout$y
	)

	plot_subnetwork(filtered_g, master_coords, "Global Network", output)
})

# -----------------------------------------------------------------------------
# 6. Phenotype-induced subnetworks
# -----------------------------------------------------------------------------

write("# Phenotype-induced subnetworks:", file="")

# Construct phenotype subnetworks and remove within-phenotype isolates.
# Matched random node sets are sampled from the full thresholded graph.
g_exp <- g

g1 <- make_phenotype_graph(
  g_exp,
  phenotype_nodes$Insomnia
)
g2 <- make_phenotype_graph(
  g_exp,
  phenotype_nodes$Schizophrenia
)
g3 <- make_phenotype_graph(
  g_exp,
  phenotype_nodes$Neuroticism
)

phenotype_networks <- list(
  Schizophrenia = g2,
  Insomnia = g1,
  Neuroticism = g3
)

# Audit exactly which phenotype-associated nodes survive the within-phenotype
# edge restriction and isolate deletion.
phenotype_subnetwork_node_audit <- bind_rows(lapply(
  names(phenotype_networks),
  function(pheno) {
    all_nodes <- phenotype_nodes[[pheno]]
    used_nodes <- V(phenotype_networks[[pheno]])$name

    data.frame(
      phenotype = pheno,
      phenotype_nodes_before_isolate_removal = length(all_nodes),
      nodes_used_in_subnetwork = length(used_nodes),
      nodes_removed_as_isolates = length(setdiff(all_nodes, used_nodes)),
      stringsAsFactors = FALSE
    )
  }
))

phenotype_null_results <- list()
for (pheno in c("Schizophrenia", "Insomnia", "Neuroticism")) {
  sub_nodes <- V(phenotype_networks[[pheno]])$name

  phenotype_null_results[[pheno]] <- integrated_subnetwork_significance(
    g_exp,
    sub_nodes = sub_nodes,
    n_sim = N_SUBNETWORK_NULL,
    indicator_attr = "indicator",
    exact_indicator_levels = FALSE,
    include_entropy = TRUE,
    alternative = "two.sided",
    p_adjust_method = "BH"
  )
}

subnetwork_results <- bind_rows(lapply(
  names(phenotype_null_results),
  function(pheno) {
    phenotype_null_results[[pheno]]$results %>%
      mutate(phenotype = pheno, .before = 1)
  }
))

print(subnetwork_results %>% mutate(across(where(is.numeric), ~ round(., 4))))

# TableS7
write.csv(
  subnetwork_results %>% mutate(across(where(is.numeric), ~ round(., 4))),
  file.path(OUTPUT_DIR, "phenotype_subnetwork_null_comparisons_TableS7.csv"),
  row.names = FALSE
)

ComparePhenotypeCentrality(g, network_cnes)
