# Multi-layered characterization of ~700,000 conserved noncoding elements within the human genome

This repository contains the core scripts used to identify and characterize conserved noncoding elements (CNEs) as described in:

“Multi-layered characterization of ~700,000 conserved noncoding elements within the human genome”

Silvia Fibi-Smetana, Fernando Fernández-Mendoza, Leila Taher

Preprint available at: https://doi.org/10.64898/2026.06.16.732547

If you use adaptations of any of these scripts or pipelines as part of your work, please cite our article. 

Please read the Materials and Methods section of the manuscript before proceeding.


## Introduction


The repository contains custom analytical scripts, lightweight example datasets, and clear instructions for reproducing the analyses. <span style="color: red;"><strong> The raw input data are not owned by us and therefore cannot be redistributed or re-hosted in this GitHub repository. </strong></span> These data are hosted by the respective public databases (e.g., UCSC and ENCODE) and are subject to their applicable data access, use, and redistribution terms.

<span style="color: red;"><strong> Instead, the repository documentation provides direct links to the original sources and specifies which files are required for each analysis. Users must download the raw input data directly from these external databases, in accordance with their respective terms of use, before running the analyses. </strong></span> This approach ensures that the repository contains only materials that we are permitted to distribute while still providing all information necessary to reproduce the analyses.


## Requirements

The scripts were developed for Linux environments, written and tested with the following tools:

Name                    | Version                 |
------------------------|-------------------------|
bigBedToBed       |                      v357|
BEDTools      |                      v2.30.0|
Python|			v3.12.7
pyBigWig python package|	v0.3.23|
samtools|			v1.19.2|
MafFilter|			v1.3.1|
Perl|				v5.30.0|
deepTools|			v3.5.4|
R|				v4.4.1|
R-package ape|		v5.8-1|
R-package clusterProfiler|	v4.14.6|
R-package vegan|		v2.7-2|
R-package igraph|		v2.2.1|
R-package dplyr|		v1.2.1|
R-package tidyr|		v1.3.2|
R-package ggplot2|		v4.0.2|
R-package gridExtra|		v2.3|
R-package cowplot|		v1.2.0|
R-package RColorBrewer|	v1.1-3|
R-package fs|			v1.6.6|
R-package data.table|		v1.18.0|
R-package stringr|		v1.6.0|
R-package scales|		v1.4.0|
R-package colorspace|	v2.1-2|
R-package RRphylo|		v3.0.2|
R-package Matrix|		v1.7-1|
R-package pheatmap|		v1.0.13|
R-package aplot|		v0.2.9|
R-package ggtree|		v3.14.0|
R-package phytools|		v2.5-2|
R-package tibble|		v3.3.1|
R-package purrr|		v1.2.1|
R-package ggrepel|		v0.9.6|

The functional characterization of the CNEs **requires downloading and pre-processing data from the ENCODE data repository**. The script `download_from_ENCODE.sh` is provided for your convenience. Before starting with the analysis, run from the project folder:

- `./data/download_from_ENCODE.sh`

This script will download the ENCODE files listed in:

- `./data/ENCODE_SCREEN_list_DNase.txt`
- `./data/ENCODE_SCREEN_list_H3K27ac.txt`

Note that “#” is used to comment out datasets.

As the downloaded bam files are rather large, they are directly converted into bigWig format and the bam files are deleted within the script.

Note that this step is **time and space intensive**.

Additionally, as described in the corresponding sections, SCREEN and GWAS datasets have to be downloaded before the functional characterization and network analysis, respectively.

## Repository Structure

The repository is organized into the following folders:

- `data`: Contains input data from external sources required for all calculations.
- `scripts`: Contains the scripts necessary to reproduce the core analyses in the paper.
- `output`: Contains results files.

## Usage:

Run the main.sh script from the main project folder: 

`./scripts/main.sh`

This script depends on the following configuration file, which contains variable names used across all scripts):

`./scripts/arguments_main_script_CNEcharacterization.cfg`

It calls the following bash scripts in the `scripts` folder:

1. `check_extract_CNEs_v2.sh`: Identifies conserved elements starting from PhastCons.
2. `annotate_phastconselements_v2.sh`: Annotates conserved elements as coding and noncoding.
3. `extract_orthologs_v3.sh`: Extracts CNEs orthologs from the multiple sequence alignment associated with the PhastCons elements and performs the ancestor estimation.
4. `functional_analysis_v2.sh`: Functionally characterizes the CNEs.
5. `network_analysis.sh`: Performs the network analysis.


### Conserved Element Identification

**Input:**

The original analysis uses the PhastCons elements from the UCSC 470-way mammalian alignment (hg38), which can be downloaded in BigBed format from https://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons470way/hg38.phastConsElements470way.bb (last accessed on September 2026)

- A reduced dataset is provided in the `data` folder (`hg38.phastConsElements470way_demo.bb`).

**Description:**

This script `check_extract_CNEs_v2.sh` depends on the following scripts:

- `merge_phastcons_codition_v2.pl`: Merges neighboring CNEs with less than 50 bp when the resulting element has more than 50% PhastCons content.

**Output:**

- `output/00_cneoutput/hg38.phastConsElements470way_demo_incl_selected_CElength.txt`: All merged elements longer than 50 bp.

### Conserved Element Annotation

**Input:**

The original analysis uses https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.chr.annotation.gff3.gz (last accessed on September 2026)
- A reduced dataset is provided in the `data` folder (`gencode.v45.annotation_demo.gff3.gz`).

**Description:**

The script `annotate_phastconselements_v2.sh` depends on the following scripts:

- `phastovmult_basic_annotation.R`: Separates elements into coding and noncoding or noncoding RNAs and other noncoding elements.
- `phastovmult_cne_annotation.R`: Defines intronic and intergenic CNEs.

**Output:**

- `output/01_annotationoutput/noncoding`: Text file containing the coordinates of CNEs.

### Identification of CNE Orthologs and CNE Ancestral Reconstruction

**Input:**

The original analysis uses:

- Multiple Sequence Alignment (MSA) files for all standard human chromosomes, available at alignment files https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz470way/maf/ (last accessed on September 2026). 
    - A reduced dataset is provided in the `data` folder (`chr1_demo.maf.gz`)

- Chromosome/Scaffold Information: Obtained for involved assemblies via:
    - UCSC chromInfo tables (e.g., http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chromInfo.txt.gz, last accessed on September 2026).
    - NCBI database using GCA accession IDs (e.g., https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_006542625.1/download?include_annotation_type=SEQUENCE_REPORT, last accessed on September 2026).
     - running faidx (v1.19.2) from SAMtools on assemblies downloaded from the DNA ZOO Project (e.g. https://dnazoo.s3.wasabisys.com/Macaca_fuscata/Macaca_fuscata_HiC.fasta.gz, last accessed on September 2026).
    - For your convenience, these data are provided in the `data/genomes_chr_length` folder.

- A configuration/options file that specifies the input, output, and sequence of filters for mafFilter.
    - Provided in the `data` folder (`maffilter_option_file`)
- A phylogenetic tree model file used alongside the multiple sequence alignment to compute conservation scores and identify PhastCons elements.
    - Provided in the `data` folder (`hg38.phyloP470way.mod`)
- A tab-separated file containing a numbered list of the assemblies used for the ortholog extraction. Assembly names should match those in the phylogenetic tree and the MSA files. Column 1: assembly number, starting from zero, Column 2: assembly name.
    - Provided in the `data` folder (`list_of_genomes.txt`)
- A list of the assemblies that are used for the ancestral reconstruction.
    - Provided in the `data` folder (`N50_filtered_species`)
- A tab-separated table specifying ancestral clades and their corresponding node numbers in the phylogenetic tree (Note: Internal node numbers must follow R's ape package indexing convention, starting at $N + 1$ where $N$ is the number of tree tips). Column 1: clade name, Column 2: node number.
    - Provided in the `data` folder (`group_nodes.txt`)


**Description:**

The script `extract_orthologs_v3.sh` depends on the following scripts:

- `my_maftobed_v2.pl`: Transforms the alignment block information into bed format.
- `extract_orthologs.pl`: Puts the ortholog information of each assembly in a separate line in the bed file.
- `change_neg_strand_positions_to_pos_strand.pl`: Calculates the positive strand coordinates for all orthologs that were on the negative strand in their assembly.
- `merge_to_phastcon_v2.pl`: Merges orthologs with consecutive information from different alignment blocks into one ortholog.
- `fit_orthologs.pl`: Corrects the ortholog coordinates from the blocks for the overlap with the CNEs.
- `lineup_orthologs_v2.pl`: Merges orthologs with less than 5bp distance in between. 
- `orth_filtering_bestcoverageperphastcon_v2.pl`: Identifies the ortholog with the largest alignment span among multiple candidate orthologs.
- `add_only_set_numbers_v2.pl`: Adds an identifier to each CNE and its orthologs.
- `demo_tree_estimation_v2.R`: Performs ancestral reconstruction.
- `transpose.awk`: Transposes a table.
- `extract_est_vector_v2.R`: Get the set numbers for the predicted sets in the different common ancestors.
- `extract_est_sets_v2.pl`: Extracts the CNE lines for the extracted set numbers.
- `combine_CNE_table_v2.R`: Combines the CNEs from different clades in a single file.

**Output:**

- `output/02_mafoutput/all_estimated_CNEs`: CNEs used for ancestral reconstruction, in bed format.
- `output/02_mafoutput/Mammals_estimated_CNEs`: Mammalian CNEs, in GTF format .
- `output/02_mafoutput/Primates_estimated_CNEs`: Primate CNEs, in GTF format .
- `output/02_mafoutput/summary_CNEs_in_clades.txt`: Tab-separated file with CNEs in rows and the clades in which they are present in the columns. Cells contain 1 if the element is conserved in the clade and NA if it is not. This file contains the data reported in **Supplementary Table 1** of the publication.


### Functional Analysis

**Input:**

- A file listing SCREEN cell type-agnostic datasets of promoters and enhancers, provided as `./data/summary_SCREEN_files.txt`, and the corresponding datasets. To download the datasets, access https://screen.encodeproject.org/ (last accessed on September 2026) from your favorite Web browser, go to “Downloads”, click on “Download human candidate promoters (hg38)” and “Download human candidate enhancers (hg38)”, and store the files in:
     `./downloads/ENCODE/SCREEN/Combined`
- A file listing SCREEN cell type-specific datasets  of candidate regulatory elements, provided as `./data/single_SCREEN_files.txt`, and the corresponding datasets. To download the datasets, access https://screen.encodeproject.org/ (last accessed on September 2026), go to “Downloads”, click on “Detailed Elements”, and search for the corresponding human datasets under “Download cCREs by cell type”. Store the files in:
     `./downloads/ENCODE/SCREEN/Single`
- A list of hg38 chromosome length formatted for BEDTools with two tab-separated columns. First column the chromosome, second column the chromosome length.
    - Provided in the `data/genomes_chr_length/hg38` folder (`chromInfo_for_betools.txt`)
- A tab separated table of ENCODE-IDs for DNase analysis with the tissue in the first column and the ENCODE-ID in the second column.
    - Provided in the `data` folder (`ENCODE_SCREEN_list_DNase.txt`)
- A tab-separated table of ENCODE-IDs for H3K27ac analysis with the tissue in the first column, first alignment file ENCODE-ID in the second column, second alignment file ENCODE-ID in the third column and control file ENCODE-ID in the fourth column.
    - Provided in the `data` folder (`ENCODE_SCREEN_list_H3K27ac.txt`)
- Transcription start sites from BioMart:
    - Provided in the `data` folder (`mart_export_TSS_Feb2025_relevant.bed`)


**Description:**

This script `functional_analysis_v2.sh` depends on the following scripts:

- `plot_sum_fisher_v2.R`: Plots the results of the summary SCREEN files bedtools fisher analysis.
- `plot_single_fisher_v2.R`: Plots the results of the single SCREEN files bedtools fisher analysis.
- `getbigwigsignalforbed_v1.py`: Takes a coordinates file and a bigwig file and gives the mean signal for each set of coordinates.
- `select_clade.R`: Splits CNE files with mean signals into mammalian and primate files.
- `DNase_mean_to_fc_v2.R`: Calculates fold-changes from the mean signals.
- `add_category_fc_v2.R`: Adds five categories based on the fold-change.
- `EpigeneticCategory_analysis_v2.R`: Calculates the number of CNEs per DNase/H3K27ac category in each tissue and defines accessible/acetylated CNEs.
- `Enrichment_v2.R`: Calculates enrichment of accessible/acetylated mammalian and primate CNEs compared to random regions.
- `correct_means_v2.R`: Computes fold-changes from the H3K27ac signals.
- `combine_categories_v2.R`: Identifies CNEs with categories 4 and 5 for accessibility and acetylation as active in the corresponding tissue.
- `Consensus_analysis_v3.R`: Counts the number of active CNEs per category and tissue.
- `Saturation_analysis_v2.R`: Performs a saturation analysis of active CNEs across tissues.
- `get_CNE_w-wo_TSS_promoter_v2.R`: Decides which CNEs are potential promoters and which not.

**Output:**

- `output/03_functionaloutput/random_estimated_CNEs`: Tab-separated file with the random regions corresponding to all_estimated_CNEs
- `output/03_functional_output/fisher_sum_summary.txt_SumFisher_barchart.pdf: Stacked barchart of clade versus Fisher´s exact test odds-ratio. Color represents the cell type-agnostic SCREEN file.
- `output/03_functional_output/fisher_sum_summary.txt_SumFisher_dotplot.pdf`: Dotplot with the number of cell type-agnostic SCREEN elements versus the number of CNEs that are overlapping with them. Color represents the clade and size the Fisher´s exact test odds-ratio.
- `output/03_functional_output/fisher_single_summary.txt_barchart.pdf`: Stacked barchart of clade versus Fisher´s exact test odds-ratio. Color represents the cell type-specific SCREEN file.
- `output/03_functional_output/fisher_single_summary.txt_dotplot.pdf`: Dotplot with the number of cell type-specific SCREEN elements versus the number of CNEs that are overlapping with them. Color represents the clade and size the Fisher´s exact test odds-ratio.
- `output/03_functionaloutput/Mammals_estimated_CNEs_small_prom_overlap`: Tab-separated table with CNE chr, start, end in first three columns and ENSG-ID in the fourth column when there was an overlap for the CNE and “.” if not.
- `output/03_functionaloutput/Primates_estimated_CNEs_small_prom_overlap`: tab-separated table where the first three columns specify the chromosome, start position, and end position of each CNE. The fourth column lists the Ensembl gene ID(s) for any overlapping TSS(s), or a period (.) if no overlap is found.
- In the `output/03_functionaloutput/DNase` folder:
    - `all_estimated_CNEs__Signal.tsv`: Mean DNase fold changes for the CNEs in the rows across the tissues in the columns.
    - `all_estimated_CNEs__Categories.tsv`: DNase categories for the CNEs in the rows across the tissues in the columns. This file contains the data reported in **Supplementary Table 2** of the publication.
    - `all_estimated_CNEs__PercentInCategories.pdf`: Barplot with fractions of mammalian and primate CNEs in the five DNase categories across tissues. This plot resembles **Figure 2d** of the publication.
    - `all_estimated_CNEs__NrOfTissuesShowingActivity.pdf`: Histogram with the relative frequency of CNEs being accessible across a certain number of tissues. This plot resembles **Figure 2e** of the publication.
    - `all_estimated_CNEs__enrichment.pdf`: Enrichment of accessible mammalian and primate CNEs compared to random genomic regions. This plot resembles **Figure 2b** of the publication.
- In the `output/03_functionaloutput/H3K27ac` folder file names are analogous to those in `output/03_functionaloutput/DNase`, but refer to H3K27ac signals
- In the `output/03_functionaloutput/consensus_DNase_H3K27ac` folder:
    - `all_estimated_CNEs__Categories.tsv_consensusactivelist`: CNEs in rows across tissues in columns. Entry is TRUE if the CNE is active in this tissue and FALSE if not. This file contains the data reported in **Supplementary Table 4** of the publication.
    - `all_estimated_CNEs__Categories.tsv_consensusactivelist_out.txt_Complexity.pdf`: Scatter plot of the number of CNEs versus the number of tissues with activity. This plot resembles **Figure 2k** of the publication.
    - `all_estimated_CNEs__Categories.tsv_consensusactivelist_saturation.txt.pdf`: Curve of the cumulative number of active CNEs versus the number of sampled tissues. This plot resembles **Figure 2j** of the publication.


### Network Analysis

**Input:**

The network analysis should be run on the full dataset (or at least, a larger portion of it) to be meaningful (and possible: it cannot be performed if no active CNEs overlap with relevant GWAS variants). The input files necessary for the analysis are generated by the aforementioned scripts:
    - `output/03_functionaloutput/Mammals_estimated_CNEs_small_prom_overlap` and `output/03_functionaloutput/Primates_estimated_CNEs_small_prom_overlap` (see “Functional Analysis”).
   - `output/03_functionaloutput/consensus_DNase_H3K27ac/all_estimated_CNEs__Categories.tsv_consensusactivelist` (see “Functional Analysis”). 
   - `output/02_mafoutput/summary_CNEs_in_clades.txt` (see “Identification of CNE Orthologs and CNE Ancestral Reconstruction”)

- GWAS-variants for certain diseases from the GWAS Catalog. Search for “insomnia” and click “Download associations”. Save the variants in a file in the `data` folder called `Insomnia_gwas-association-downloaded_2026-02-05-EFO_0004698.tsv`. Then search for “neuroticism measurement” and save the variants in a file called `Neuroticism_measurement_gwas-association-downloaded_2026-02-05-EFO_0007660.tsv`. Search for “neuroticism”, check “Include child trait data”, download and save the variants in a file called `Neuroticism_gwas-association-downloaded_2026-02-05-EFO_0004257-withChildTraits.tsv`. Repeat for “schizophrenia”; save the data as  `Schizophrenia_gwas-association-downloaded_2026-02-05-MONDO_0005090-withChildTraits.tsv`.

**Description:**

This script `network_analysis.sh` depends on the following script:
- `network_analyses.R`: Performs all network analyses described in the manuscript.


**Output:**

- `output/04_networkoutput/pruned_network_metrics_TableS6.csv`:  Topological properties of the pruned CNE co-activity network. Comma-separated table with the metric in the first column and the values for mean, SD and median in columns two to four. This file contains the data reported in **Supplementary Table 6** of the publication.
- `output/04_networkoutput/phenotype_subnetwork_null_comparisons_TableS7.csv`: . Topological properties of the schizophrenia, insomnia and neuroticism induced subnetworks. Comma-separated table with phenotype in the first column, the statistical descriptor in the second column and the values (observed, null_mean, null_median, null_sd, p-value, p_adj, z) in columns three to nine. This file contains the data reported in **Supplementary Table 7** of the publication.
- `output/04_networkoutput/network_analysis.R_NetworkGraph.pdf`: Visual representation of the global network. This plot resembles **Figure 4c** of the publication.
- `output/04_networkoutput/network_analysis.R_NetworkPrunning.pdf`: Visual representation of the topological changes in the CNE co-activity network upon recursive weight pruning. This plot resembles **Figure S11** of the publication.
- `output/04_networkoutput/network_analysis.R_NetworkCNEConservationByTrait.pdf`: Stacked barplot of the number of CNEs for each trait belonging to promoters and enhancers and mammalian and primate clade. This plot resembles **Figure 4b** of the publication.
- `output/04_networkoutput/network_analysis.R_CentralityVsConservation`: Centrality profile distribution of enhancer-annotated CNEs by evolutionary conservation. This plot resembles **Figure 4e** of the publication.
- `output/04_networkoutput/percolation_analysis.pdf`: Visual representation of the network percolation.


## Reproducibility

Because the analyses rely on:

* Large multiple sequence alignments
* ENCODE datasets
* Dataset-specific preprocessing steps

The scripts provided here capture the **core analytical logic** used in the study. Full reproduction of the results requires access to the original data and following the workflow described in the manuscript.
