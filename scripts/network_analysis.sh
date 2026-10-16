#!/bin/bash

# Input Data

[[ -f "$ALLESTSETS" ]] || error_exit "$LINENO" "Cannot find $ALLESTSETS!"

[[ -f "$SCHIZOTSV" ]] || error_exit "$LINENO" "Cannot find $SCHIZOTSV!"
[[ -f "$INSOTSV" ]] || error_exit "$LINENO" "Cannot find $INSOTSV!"
[[ -f "$NEUROTSV" ]] || error_exit "$LINENO" "Cannot find $NEUROTSV!"
[[ -f "$NEUROMTSV" ]] || error_exit "$LINENO" "Cannot find $NEUROMTSV!"

# CNE Activity

CNE_ACTIVITY="${CONSDIR}/"$(basename "$ALLESTSETS")"__Categories.tsv_consensusactivelist"
[[ -f "$CNE_ACTIVITY" ]] || error_exit "$LINENO" "Cannot find $CNE_ACTIVITY!"

# CNE Conservation

[[ -f "$SUM" ]] || error_exit "$LINENO" "Cannot find $SUM!"

## CNEs overlapping with Promoters

PRIMATE_PROMOTER_FILE="${FUNCTIONALOUTPUTDIR}/Primates_estimated_CNEs_small_prom_overlap"
[[ -f "$PRIMATE_PROMOTER_FILE" ]] || error_exit "$LINENO" "Cannot find $PRIMATE_PROMOTER_FILE!"
MAMMAL_PROMOTER_FILE="${FUNCTIONALOUTPUTDIR}/Mammals_estimated_CNEs_small_prom_overlap"
[[ -f "$MAMMAL_PROMOTER_FILE" ]] || error_exit "$LINENO" "Cannot find $MAMMAL_PROMOTER_FILE!"

# Output

[[ -d "$NETWORKOUTPUTDIR" ]] || error_exit "$LINENO" "Cannot find $NETWORKOUTPUTDIR!"


######################################################################################
# GWAS DISEASES
######################################################################################

###Overlap CNEs and variants in selected disaeses
echo "# Overlap CNEs and variants in selected disaeses"

SCHIZOBED="${DEMOFILEDIR}/Schizophrenia_gwas.bed"
INSOBED="${DEMOFILEDIR}/Insomnia_gwas.bed"
NEUROBED="${DEMOFILEDIR}/Neuroticism_gwas.bed"
NEUROMBED="${DEMOFILEDIR}/Neuroticismmeasurement_gwas.bed"

tmp1=$(mktemp -u -p .)

# Schizophrenia
input="$SCHIZOTSV"
ALLESTSETSSCHIZO="${NETWORKOUTPUTDIR}/"$(basename "$ALLESTSETS")"_Schizophrenia_new_intersected_filtered"
echo "# $input"
cat "$input" | tail -n +2 | awk -F'\t' '{print $12"\t"$13"\t"$13+1}' | sed -e 's/^/chr/' | sort -k1,1V -k2,2n | awk -F'\t' '{if($3!=1) print $0}' > "$tmp1"
head -n3 "$tmp1"
wc -l "$tmp1"
bedtools intersect -a "$ALLESTSETS" -b "$tmp1" -wao | awk -F"\t" '{if($4!=".") print $1"\t"$2"\t"$3}' | uniq > "$ALLESTSETSSCHIZO"
head -n3 "$ALLESTSETSSCHIZO"
wc -l "$ALLESTSETSSCHIZO"
rm "$tmp1"

# Insomnia
input="$INSOTSV"
ALLESTSETSINSO="${NETWORKOUTPUTDIR}/"$(basename "$ALLESTSETS")"_Insomnia_new_intersected_filtered"
echo "# $input"
cat "$input" | tail -n +2 | awk -F'\t' '{print $12"\t"$13"\t"$13+1}' | sed -e 's/^/chr/' | sort -k1,1V -k2,2n | awk -F'\t' '{if($3!=1) print $0}' > "$tmp1"
head -n3 "$tmp1"
wc -l "$tmp1"
bedtools intersect -a "$ALLESTSETS" -b "$tmp1" -wao | awk -F"\t" '{if($4!=".") print $1"\t"$2"\t"$3}' | uniq > "$ALLESTSETSINSO"
head -n3 "$ALLESTSETSINSO"
wc -l "$ALLESTSETSINSO"
rm "$tmp1"

# Neuroticism

input="$NEUROTSV"
ALLESTSETSNEURO="${NETWORKOUTPUTDIR}/"$(basename "$ALLESTSETS")"_Neuroticism_new_intersected_filtered"
echo "# $input"
cat "$input" | tail -n +2 | awk -F'\t' '{print $12"\t"$13"\t"$13+1}' | sed -e 's/^/chr/' | sort -k1,1V -k2,2n | awk -F'\t' '{if($3!=1) print $0}' > "$tmp1"
head -n3 "$tmp1"
wc -l "$tmp1"
bedtools intersect -a "$ALLESTSETS" -b "$tmp1" -wao | awk -F"\t" '{if($4!=".") print $1"\t"$2"\t"$3}' | uniq > "$ALLESTSETSNEURO"
head -n3 "$ALLESTSETSNEURO"
wc -l "$ALLESTSETSNEURO"
rm "$tmp1"

# Neuroticism Measurement

ALLESTSETSNEUROM="${NETWORKOUTPUTDIR}/"$(basename "$ALLESTSETS")"_Neuroticismmeasurement_new_intersected_filtered"
echo "# $input"
cat "$input" | tail -n +2 | awk -F'\t' '{print $12"\t"$13"\t"$13+1}' | sed -e 's/^/chr/' | sort -k1,1V -k2,2n | awk -F'\t' '{if($3!=1) print $0}' > "$tmp1"
head -n3 "$tmp1"
wc -l "$tmp1"
bedtools intersect -a "$ALLESTSETS" -b "$tmp1" -wao | awk -F"\t" '{if($4!=".") print $1"\t"$2"\t"$3}' | uniq > "$ALLESTSETSNEUROM"
head -n3 "$ALLESTSETSNEUROM"
wc -l "$ALLESTSETSNEUROM"
rm "$tmp1"

######################################################################################
# NETWORK ANALYSIS
######################################################################################

echo "# Network analysis"

[[ -f "$ALLESTSETSINSO" ]] || error_exit "$LINENO" "Cannot find $ALLESTSETSINSO!"
[[ -f "$ALLESTSETSSCHIZO" ]] || error_exit "$LINENO" "Cannot find $ALLESTSETSSCHIZO!"
[[ -f "$ALLESTSETSNEURO" ]] || error_exit "$LINENO" "Cannot find $ALLESTSETSNEURO!"
[[ -f "$ALLESTSETSNEUROM" ]] || error_exit "$LINENO" "Cannot find $ALLESTSETSNEUROM!"

read -p "The network analysis will be performed on the full dataset for demonstration. Comment the variables below if you would prefer to run it otherwise (Note: This is not guaranteed to work and may fail if there is not enough data.). Confirm [y/N]: " confirm

case "$confirm" in
    [yY][eE][sS]|[yY])
        echo "Proceeding..."
        # Add your script logic here
        ;;
    *)
        echo "Aborted."
        exit 1
        ;;
esac

# Uncomment to run the script on the full dataset
FULLDATASET="./network_full_dataset"

CNE_ACTIVITY="$FULLDATASET/all_estimated_CNEs__Categories.tsv_consensusactivelist"
[[ -f "$CNE_ACTIVITY" ]] || error_exit "$LINENO" "Cannot find $CNE_ACTIVITY!"
SUM="$FULLDATASET/summary_CNEs_in_clades.txt"
[[ -f "$SUM" ]] || error_exit "$LINENO" "Cannot find $SUM!"

PRIMATE_PROMOTER_FILE="$FULLDATASET/Primates_estimated_CNEs_small_prom_overlap"
[[ -f "$PRIMATE_PROMOTER_FILE" ]] || error_exit "$LINENO" "Cannot find $PRIMATE_PROMOTER_FILE!"
MAMMAL_PROMOTER_FILE="$FULLDATASET/Mammals_estimated_CNEs_small_prom_overlap"
[[ -f "$MAMMAL_PROMOTER_FILE" ]] || error_exit "$LINENO" "Cannot find $MAMMAL_PROMOTER_FILE!"

ALLESTSETSINSO="$FULLDATASET/all_estimated_CNEs_Insomnia_new_intersected_filtered"
[[ -f "$ALLESTSETSINSO" ]] || error_exit "$LINENO" "Cannot find $ALLESTSETSINSO!"
ALLESTSETSSCHIZO="$FULLDATASET/all_estimated_CNEs_Schizophrenia_new_intersected_filtered"
[[ -f "$ALLESTSETSSCHIZO" ]] || error_exit "$LINENO" "Cannot find $ALLESTSETSSCHIZO!"
ALLESTSETSNEURO="$FULLDATASET/all_estimated_CNEs_Neuroticism_new_intersected_filtered"
[[ -f "$ALLESTSETSNEURO" ]] || error_exit "$LINENO" "Cannot find $ALLESTSETSNEURO!"
ALLESTSETSNEUROM="$FULLDATASET/all_estimated_CNEs_Neuroticismmeasurement_new_intersected_filtered"
[[ -f "$ALLESTSETSNEUROM" ]] || error_exit "$LINENO" "Cannot find $ALLESTSETSNEUROM!"

R --slave --no-save --no-restore --no-environ --silent --args "$CNE_ACTIVITY" "$SUM" "$PRIMATE_PROMOTER_FILE" "$MAMMAL_PROMOTER_FILE" "$ALLESTSETSINSO" "$ALLESTSETSSCHIZO" "$ALLESTSETSNEURO" "$ALLESTSETSNEUROM" "$NETWORKOUTPUTDIR" < "${SCRIPTS}/network_analyses.R"
