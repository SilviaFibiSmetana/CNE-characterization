#!/bin/bash

set -Eeuo pipefail

######################################################
##This script calls all the other necessary scripts
##for the complete analysis of this project
##
######################################################

trap cleanup SIGINT SIGTERM ERR EXIT;

error_exit(){
    line=$1;
    shift 1;
    echo "ERROR: non zero return code from line: $line - $@" >&2;
    exit 1;
}

cleanup(){
    trap - SIGINT SIGTERM ERR EXIT
}

# Load the configuration file

echo "# Loading the configuration file..."

WORKDIR="$PWD/"

echo "# Working directory: $WORKDIR"

source "${WORKDIR}scripts/arguments_main_script_CNEcharacterization.cfg"

# Check input files

[[ -f "$PHASTCONSELEMENTSFILEBB" ]] || error_exit "$LINENO" "Cannot find $PHASTCONSELEMENTSFILEBB!";
[[ -f "$GENOMEGFF" ]] || error_exit "$LINENO" "Cannot find $GENOMEGFF!";

echo "# Using PhastCons elements in $PHASTCONSELEMENTSFILEBB"
echo "# Using genome annotation in $GENOMEGFF"
echo "# Main output directory: $OUTPUTDIR"

[[ -d "$OUTPUTDIR" ]] || mkdir "$OUTPUTDIR";

# Create output directories for the different parts of the script

[[ -d "$CNEOUTPUTDIR" ]] || mkdir "$CNEOUTPUTDIR"
[[ -d "$ANNOTATIONOUTPUTDIR" ]] || mkdir "$ANNOTATIONOUTPUTDIR"
[[ -d "$MAFOUTPUTDIR" ]] || mkdir "$MAFOUTPUTDIR"
[[ -d "$FUNCTIONALOUTPUTDIR" ]] || mkdir "$FUNCTIONALOUTPUTDIR"
[[ -d "$NETWORKOUTPUTDIR" ]] || mkdir "$NETWORKOUTPUTDIR"

# Main Analysis

echo "
######################################################
#
# 00 Identification of Conserved Elements from PhastCons
#
######################################################
"

. "${SCRIPTS}/check_extract_CNEs_v2.sh"


echo "
######################################################
#
# 01 Conserved Element Annotation
#
######################################################
"

. "${SCRIPTS}/annotate_phastconselements_v2.sh"


echo "
######################################################
##
## 02 Ortholog Extraction and Ancestral Reconstruction
##
######################################################
"

. "${SCRIPTS}/extract_orthologs_v3.sh"


echo "
######################################################
##
## 03 CNE Functional Characterization
##
######################################################
"

. "${SCRIPTS}/functional_analysis_v2.sh"


echo "
######################################################
##
## 04 CNE GWAS and Network Analysis
##
######################################################
"

. "${SCRIPTS}/network_analysis.sh"

echo "Done!"
