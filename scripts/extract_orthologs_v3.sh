#!/bin/bash

[[ -d "$CNEOUTPUTDIR" ]] || error_exit "$LINENO" "Cannot find $CNEOUTPUTDIR!"
[[ -f "$ELEMENTSNONCODING" ]] || error_exit "$LINENO" "Cannot find $ELEMENTSNONCODING!"
[[ -f "$TREE" ]] || error_exit "$LINENO" "Cannot find $TREE!"
[[ -f "$GENOMEINDEX" ]] || error_exit "$LINENO" "Cannot find $GENOMEINDEX!"
[[ -f "$MAFFILTEROPTIONFILE" ]] || error_exit "$LINENO" "Cannot find $MAFFILTEROPTIONFILE!"
[[ -f "$FILTEREDSPECIES" ]] || error_exit "$LINENO" "Cannot find $FILTEREDSPECIES!"
[[ -f "$CLADENODES" ]] || error_exit "$LINENO" "Cannot find $FILTEREDSPECIES!"

if [[ -z "${ORTHOLOGDIST:-}" ]]; then
  error_exit "$LINENO" "$ORTHOLOGDIST is undefined or empty."
fi

if ! declare -p CladesOfInterest 2>/dev/null | grep -q 'declare -a'; then
	error_exit "$LINENO" "CladesOfInterest does not exist or is not an array."
fi

DISTDIR="${MAFOUTPUTDIR}/orthologs_linedup_${ORTHOLOGDIST}_as_dist"
[[ -d "$DISTDIR" ]] || mkdir "$DISTDIR"

tmp1=$(mktemp -u -p .)
tmp2=$(mktemp -u -p .)
tmp3=$(mktemp -u -p .)
tmp4=$(mktemp -u -p .)


# PREPARE THE LIST OF ELEMENTS TO EXTRACT THEM FROM THE MAF FILES

FEATURES="${MAFOUTPUTDIR}/noncoding_for_maffileextraction.gtf"

cat "$ELEMENTSNONCODING" | awk '{print $1"\t""\.""\t""\.""\t"$2"\t"$3"\t""\.""\t""\+""\t""\.""\t""gene_id \"test\";"}' > "$FEATURES"
head -n 5 "$FEATURES"
wc -l "$FEATURES"

# EXTRACT ALL THE CHROMOSOME NAMES FOR THE TARGET SPECIES

CHRFILE="${MAFOUTPUTDIR}/chrs.txt"
cat "$FEATURES" | awk '{print $1}' | uniq > "$CHRFILE"
head -n 5 "$CHRFILE"
wc -l "$CHRFILE"

# EXTRACT THE FEATURES FOR EACH CHROMOSOME 

while read -r CHR; do
	FEATURESELFILE="${MAFOUTPUTDIR}/${CHR}_features.gtf"
	awk -v chr="$CHR" '$0~chr{print}' "$FEATURES" > "$FEATURESELFILE"
	wc -l "$FEATURESELFILE"
done < "$CHRFILE"

# USE MAFFILTER TO EXTRACT THE RELEVANT BLOCKS FOR EACH CHROMOSOME ACCORDING TO THE FEATURES TO MAKE MAFFILES SMALLER

cp "$MAFFILTEROPTIONFILE" "${MAFOUTPUTDIR}/maffilter_option_file"

while read -r CHR; do
	FEATURESELFILE="${MAFOUTPUTDIR%}/${CHR}_features.gtf"
	MF_OPTION_FILE="${MAFOUTPUTDIR}/${CHR}_maffilter_option_file"
	# Clean paths because MafFilter interprets // as comments
	MF_OPTION_FILE="${MF_OPTION_FILE//\/\//\/}"
	cat "${MAFOUTPUTDIR}/maffilter_option_file" | sed "s#PATHTOCHROMOSOMENUMBER#$MAFDIR/$CHR#" | sed "s#FEATUREFILEASGTF#$FEATURESELFILE#" > "${MF_OPTION_FILE}" 
	# Clean paths because MafFilter interprets // as comments
	sed 's#//*#/#g' -i "${MF_OPTION_FILE}"
	
	echo
	cat "${MF_OPTION_FILE}"
	echo
	
	# Clean paths because MafFilter interprets // as comments
	command="maffilter param=${MF_OPTION_FILE}"
	echo "# Running $command..."
	eval "$command"
	mv "${MAFDIR}/${CHR}_extracted.maf.gz" "${MAFOUTPUTDIR}/${CHR}_extracted.maf.gz"
	mv "${MAFDIR}/${CHR}.maffilter.log" "${MAFOUTPUTDIR}/${CHR}.maffilter.log"
done < "$CHRFILE"


while read -r CHR; do
	echo "# Processing $CHR..."
	
	# Convert MAF to BED format
	MAFEXTFILE="${MAFOUTPUTDIR}/${CHR}_extracted.maf.gz"
	MAFEXTBEDFILE="${MAFOUTPUTDIR}/${CHR}_extracted.maf.bed"	
	perl "${SCRIPTS}/my_maftobed_v2.pl" "$MAFEXTFILE" "$GENOMEINDEX" "$MAFEXTBEDFILE"
	MAFEXTBEDFILEGZ="${MAFOUTPUTDIR}/${CHR}_extracted.maf.bed.gz"
	cat "$MAFEXTBEDFILE" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' | gzip -c > $MAFEXTBEDFILEGZ
	rm "$MAFEXTFILE"

	# Intersect the CNEs with the extracted maf information
	FEATURESELFILE="${MAFOUTPUTDIR}/${CHR}_features.gtf"
	FEATURESELFILEBED="${MAFOUTPUTDIR}/${CHR}_features.bed"
	awk '{print $1"\t"$4"\t"$5}' "$FEATURESELFILE" > "$FEATURESELFILEBED"
	INTERSECTED="${MAFOUTPUTDIR}/${CHR}_intersected"
	bedtools intersect -wo -a "$FEATURESELFILEBED" -b "$MAFEXTBEDFILEGZ" > "$INTERSECTED"	

	# Extract the info for each species in a single line which has a match in the maf block
	ORTHEXTR="${MAFOUTPUTDIR}${CHR}_orthextracted"
	perl "${SCRIPTS}/extract_orthologs.pl" "$INTERSECTED" > "$ORTHEXTR" 
	rm "$INTERSECTED"
	wc -l "$ORTHEXTR"

	# Remove the lines of the species that are not in the $filteredspecies list

	ORTHEXTRFILTERED="${MAFOUTPUTDIR}/${CHR}_orthextracted_filtered"
	grep -Fwf "$FILTEREDSPECIES" "$ORTHEXTR" > "$ORTHEXTRFILTERED"
	head -n3 "$ORTHEXTRFILTERED"
	rm "$ORTHEXTR"
	wc -l "$ORTHEXTRFILTERED"
	echo

	# For query sequences or alignment blocks mapped to the negative strand (such as in chain or MAF formats), local coordinates are measured backward from the end of the sequence (relative to the sequence size

	ORTHEXTRSTRANDCORR="${MAFOUTPUTDIR}/${CHR}_orthextracted_strandcorrected"
	perl "${SCRIPTS}/change_neg_strand_positions_to_pos_strand.pl" "$FILTEREDSPECIES" "$SPECIESCHRINFODIRPATH" "$ORTHEXTRFILTERED" > "$ORTHEXTRSTRANDCORR" 
	rm "$ORTHEXTRFILTERED"
	wc -l "$ORTHEXTRSTRANDCORR"

	# Sort the files with the extracted information
        ORTHEXTRSTRANDCORRSORTED="${MAFOUTPUTDIR}/${CHR}_orthextracted_strandcorrectedsorted"
        cat "$ORTHEXTRSTRANDCORR" | sort -k8,8 -k2,2n -k5,5n -k10,10 > "$ORTHEXTRSTRANDCORRSORTED"
        rm "$ORTHEXTRSTRANDCORR"
	wc -l "$ORTHEXTRSTRANDCORRSORTED"

	# Merge the lines with consecutive info into one
        ORTHMERGED="${MAFOUTPUTDIR}/${CHR}_orthologs_merged"
        perl "${SCRIPTS}/merge_to_phastcon_v2.pl" "$ORTHEXTRSTRANDCORRSORTED" "$ORTHMERGED"
	wc -l "$ORTHMERGED"

	# Obtain reference and orthologous positions for the part of the alignment overlapping the element of interest
	ORTHMERGEDOVERLAPCORR="${MAFOUTPUTDIR}/${CHR}_orthologs_merged_overlapcorr"
        perl "${SCRIPTS}/fit_orthologs.pl" "$ORTHMERGED" "$ORTHMERGEDOVERLAPCORR"
        rm "$ORTHMERGED"
	wc -l "$ORTHMERGEDOVERLAPCORR"

	# Remove the lines where the last step did not leave at least one nucleotide for the ortholog
	ORTHMERGEDOVERLAPCORRCLEAN="${MAFOUTPUTDIR}/${CHR}_orthologs_merged_overlapcorr_clean"
        cat "$ORTHMERGEDOVERLAPCORR" | awk 'NF==16 {print $0}' > "$ORTHMERGEDOVERLAPCORRCLEAN"
        rm "$ORTHMERGEDOVERLAPCORR"
	wc -l "$ORTHMERGEDOVERLAPCORRCLEAN"

	# Sort the files with the orthologs according to phastcon and ortholog
	ORTHMERGEDOVERLAPCORRCLEANSORTED="${MAFOUTPUTDIR}/${CHR}_orthologs_merged_overlapcorr_clean_sorted"
        cat "$ORTHMERGEDOVERLAPCORRCLEAN" | sort -k8,8 -k1,1 -k2,2n -k9,9 -k10,10n > "$ORTHMERGEDOVERLAPCORRCLEANSORTED"
        rm "$ORTHMERGEDOVERLAPCORRCLEAN"
	wc -l "$ORTHMERGEDOVERLAPCORRCLEANSORTED"

	# Summarize the orthologs where the ortholog distance is not bigger than a certain threshold.
	ORTHLINEDUP="${DISTDIR}/${CHR}_orthologs_linedup"
	perl "${SCRIPTS}/lineup_orthologs_v2.pl" "$ORTHMERGEDOVERLAPCORRCLEANSORTED" "$ORTHLINEDUP" "$ORTHOLOGDIST"
	rm "$ORTHMERGEDOVERLAPCORRCLEANSORTED"
	head -n3 "$ORTHLINEDUP"
	wc -l "$ORTHLINEDUP"

	# Filter the orthologs, based on statistical criteria. Among multiple putative orthologs, chose the one with the best coverage of the element in the reference genome.
	ORTHLINEDUPFILTERED1="${DISTDIR}/${CHR}_orthologs_linedup_filtered1"
        perl "${SCRIPTS}/orth_filtering_bestcoverageperphastcon_v2.pl" "$ORTHLINEDUP" "$ORTHLINEDUPFILTERED1"
        rm "$ORTHLINEDUP"
	wc -l "$ORTHLINEDUPFILTERED1"

	# Filter for coverage/length cutoff >24 by calculating the length of the CNE and multiply by the coverage to get an effective length
	cutoff=24
        ORTHLINEDUPFILTERED2B="${DISTDIR}/${CHR}_orthologs_linedup_filtered2b"
        cat "$ORTHLINEDUPFILTERED1" | awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"($3-$2+1)*$14 }' | awk -v var="$cutoff" '{ if ($15>var) print $0}' > "$ORTHLINEDUPFILTERED2B"
        rm "$ORTHLINEDUPFILTERED1"
	wc -l "$ORTHLINEDUPFILTERED2B"

	# Sort and combine the filtered ortholog files
	ORTHLINEDUPFILTERED2BSORTED="${DISTDIR}/${CHR}_orthologs_linedup_filtered2b_sorted"
	# Sort the files that the block belonging to one human phastcon is together and it has the human entry in the first row
	cat "$ORTHLINEDUPFILTERED2B" | awk 'BEGIN {OFS="\t"} {if ($7 == "hg38") print "1", $0; else print "2", $0}' | sort -k2,2V -k3,3n -k4,4n -k1,1n | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}' > "$ORTHLINEDUPFILTERED2BSORTED" #cut -f2-
	cat "$ORTHLINEDUPFILTERED2BSORTED" | awk -F'\t' '{print NF}' | sort -nu
	rm "$ORTHLINEDUPFILTERED2B"
	wc -l "$ORTHLINEDUPFILTERED2BSORTED"

done < "$CHRFILE"
echo


# MERGE THE TABLES FOR ALL CHROMOSOMES INTO ONE
echo "# Building multi-chromosomal table..."

ALLORTHOLOGS="${DISTDIR}/all_orthologs_linedup_filtered2b_sorted"
rm -f "$ALLORTHOLOGS"
for f in "${DISTDIR}"/chr*_orthologs_linedup_filtered2b_sorted; do
	echo $f
	wc -l $f
	cat "$f" >> "$ALLORTHOLOGS" 
done
wc -l "$ALLORTHOLOGS"

# ADD SET NUMBERS TO THE TABLE

ALLORTHOLOGSNUMBERED="${DISTDIR}/all_orthologs_linedup_filtered2b_sorted_numbered"
ALLORTHOLOGSSHORT="${DISTDIR}/all_orthologs_linedup_filtered2b_sorted_short"
perl "${SCRIPTS}/add_only_set_numbers_v2.pl" "$ALLORTHOLOGS" "$ALLORTHOLOGSNUMBERED"
cat "$ALLORTHOLOGSNUMBERED" | awk '{print $7"\t"$15}' > "$ALLORTHOLOGSSHORT"

# EXTRACT THE HUMAN CNES AS BED FILE

CNESBEFOREESTIMATIONS="${MAFOUTPUTDIR}/cnes_before_estimations.bed"
cat "$ALLORTHOLOGSNUMBERED" | grep "hg38" | awk '{ print $1"\t"$2"\t"$3}' | sort -k1,1V -k2,2n > "$CNESBEFOREESTIMATIONS"
head "$CNESBEFOREESTIMATIONS"
wc -l "$CNESBEFOREESTIMATIONS"
echo

# ANCESTRAL RECONSTRUCTION
echo "# Performing ancestral estimation..."

NODEESTIMATIONS="${MAFOUTPUTDIR}/all_orthologs_linedup_filtered2b_sorted_short_nodeestimations.txt"
command="R --slave --no-save --no-restore --no-environ --silent --args $TREE $MAFOUTPUTDIR $ALLORTHOLOGSSHORT $FILTEREDSPECIES < ${SCRIPTS}/demo_tree_estimation_v2.R"
echo "$command"
eval "$command"
head -n3 "$NODEESTIMATIONS"
wc -l "$NODEESTIMATIONS"
echo

# PICK JUST THE HUMAN LINE FOR THE SETS PREDICTED IN A CERTAIN CLADE 
echo "# Preparing CNE datasets..."

for clade in ${CladesOfInterest[@]}; do
	echo "# Processing $clade..."

	# Extract the node number for this clade
	NODE=$(cat "$CLADENODES" | grep "$clade" | awk '{print $2}')
	echo "$NODE"

	# Extract the set numbers from the matching sets
	CLADEESTNUMBERS="${MAFOUTPUTDIR}/${clade}_estimated_CNE_numbers.txt"
	cat "$NODEESTIMATIONS" | grep "$NODE" > "$tmp1"
	awk -f "${SCRIPTS}/transpose.awk" "$tmp1" > "$tmp2"
	command="R --slave --no-save --no-restore --no-environ --silent --args $tmp2 $clade $MAFOUTPUTDIR < ${SCRIPTS}/extract_est_vector_v2.R"
	echo "# $command"
	eval "$command"
	head -n3 "$CLADEESTNUMBERS"
	wc -l "$CLADEESTNUMBERS"
	
	# Extract the hg38 lines for the extracted set numbers
	perl "${SCRIPTS}/extract_est_sets_v2.pl" "$CLADEESTNUMBERS" "$ALLORTHOLOGSNUMBERED" "$tmp3"
	head "$tmp3"
	wc -l "$tmp3"
	cat "$tmp3" | sort -k1,1V -k2,2n > "$tmp4"
	
	# filter out CNEs shorter 2000 nucleotides
	CLADEESTSETSWITHXY="${MAFOUTPUTDIR}/${clade}_estimated_CNEsXY"
	cat "$tmp4" | awk '{if(($3-$2)<=2000) print $0}' > "$CLADEESTSETSWITHXY"
	wc -l "$CLADEESTSETSWITHXY"
	rm "$tmp1" "$tmp2" "$tmp3" "$tmp4"
	
	# Select only the CNEs on chr1-22 and remove the sex chromosomes 
	CLADEESTSETS="${MAFOUTPUTDIR}/${clade}_estimated_CNEs"
	cat "$CLADEESTSETSWITHXY" | grep -v "chrX" | grep -v "chrY" > "$CLADEESTSETS"
	wc -l "$CLADEESTSETS"
	cat "$CLADEESTSETS" | awk '{print $1}' | sort | uniq 
	echo
done

# COMBINE ALL HUMAN CNES USED FOR ANCESTRAL RECONSTRUCTION IN ONE FILE

echo "# Creating final output files..."

rm -f "$tmp1"
for clade in ${CladesOfInterest[@]}; do
	CLADEESTSETS="${MAFOUTPUTDIR}/${clade}_estimated_CNEs"
	[[ -f "$CLADEESTSETS" ]] || error_exit "$LINENO" "Cannot find $CLADEESTSETS!"
	cat "$CLADEESTSETS" | awk '{print $1"\t"$2"\t"$3}' > "$tmp1"
done
cat "$tmp1" | sort -k1,1V -k2,2n | uniq > "$ALLESTSETS"
rm "$tmp1"
head -n3 "$ALLESTSETS"
wc -l "$ALLESTSETS"


# COMBINE ALL CNES FOR WHICH ANCESTRAL RECONSTRUCTION WAS PERFORMED INTO AN OVERVIEW FILE 

MAMESTSETS="${MAFOUTPUTDIR}/Mammals_estimated_CNEs"
PRIESTSETS="${MAFOUTPUTDIR}/Primates_estimated_CNEs"
cat "$MAMESTSETS" | awk '{print $1"_"$2"_"$3}' > "$tmp1"
cat "$PRIESTSETS" | awk '{print $1"_"$2"_"$3}' > "$tmp2"
command="R --slave --no-save --no-restore --no-environ --silent --args $ALLESTSETS ${tmp1} ${tmp2} ${MAFOUTPUTDIR} < ${SCRIPTS}/combine_CNE_table_v2.R"
echo "# $command"
eval $command
rm "$tmp1" "$tmp2"

echo "# CNEs are stored in"

head -n3 "$ALLESTSETS"
wc -l "$ALLESTSETS"

head -n3 "$MAMESTSETS"
wc -l "$MAMESTSETS"

head -n3 "$PRIESTSETS"
wc -l "$PRIESTSETS"

echo "Done!"
