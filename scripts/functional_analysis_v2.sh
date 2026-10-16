#!/bin/bash

# Input Data

[[ -d "$DEMOFILEDIR" ]] || error_exit "$LINENO" "Cannot find $DEMOFILEDIR!"

[[ -f "$ALLESTSETS" ]] || error_exit "$LINENO" "Cannot find $ALLESTSETS!"
[[ -f "$HGCHRINFOBEDTOOLS" ]] || error_exit "$LINENO" "Cannot find $HGCHRINFOBEDTOOLS!"
[[ -f "$SUM" ]] || error_exit "$LINENO" "Cannot find $SUM!"

## ENCODE AND SCREEN Datasets

[[ -f "$SUMSCREENLIST" ]] || error_exit "$LINENO" "Cannot find $SUMSCREENLIST!"
[[ -f "$SINGLESCREENLIST" ]] || error_exit "$LINENO" "Cannot find $SINGLESCREENLIST!"
[[ -f "$ENCODESCREENDNaseLIST" ]] || error_exit "$LINENO" "Cannot find $ENCODESCREENDNaseLIST!"
[[ -f "$ENCODESCREENH3K27acLIST" ]] || error_exit "$LINENO" "Cannot find $ENCODESCREENH3K27acLIST!"

## TSS

[[ -f "$TSSRELEVANT" ]] || error_exit "$LINENO" "Cannot find $TSSRELEVANT!"

# Output

[[ -d "$FUNCTIONALOUTPUTDIR" ]] || error_exit "$LINENO" "Cannot find $FUNCTIONALOUTPUTDIR!"

[[ -n "${CONSDIR+x}" ]] || error_exit "$LINENO" "\$CONSDIR is not defined!"

if [[ -z "${CNE_ACTIVITY:-}" ]]; then
  error_exit "$LINENO" "$CNE_ACTIVITY is undefined or empty."
fi


# Parameters

if ! declare -p CladesOfInterest 2>/dev/null | grep -q 'declare -a'; then
	error_exit "$LINENO" "CladesOfInterest does not exist or is not an array."
fi

tmp1=$(mktemp -u -p .)
tmp2=$(mktemp -u -p .)

DATA_PREPROCESSING=true
OVERLAPS_WITH_CELLAGNOSTIC_SCREEN=true
OVERLAPS_WITH_CELLSPEC_SCREEN=true
DNASE_ANALYSIS=true
H3K27AC_ANALYSIS=true
ACTIVITY_ANALYSIS=true
OVERLAPS_WITH_PROMOTERS=true

######################################################################################
# FUNCTIONAL ANALYSIS
######################################################################################


######################################################################################
# DATA PREPARATION
######################################################################################

if [[ "$DATA_PREPROCESSING" = true ]]; then

	echo "# Data preprocessing..."
	
	## GENERATE FILES WITH ONLY THE COORDINATES OF THE CNES

	for clade in "${CladesOfInterest[@]}"; do
		echo "# Extracting CNEs for $clade..."
		cne_file="${MAFOUTPUTDIR}/${clade}_estimated_CNEs"
		[[ -f "$cne_file" ]] || error_exit "$LINENO" "Cannot find $cne_file!"
		head -n3 "$cne_file"
		wc -l "$cne_file"

		output="${MAFOUTPUTDIR}/${clade}_estimated_CNEs_small"
		cat "$cne_file" | awk '{print $1"\t"$2"\t"$3}' > "$output"
		head -n3 "$output"
		wc -l "$output"
	done	

	## SELECT A RANDOM REGION IN THE GENOME FOR EACH CNE USED FOR ANCESTRAL RECONSTRUCTION 

	wc -l "$ALLESTSETS"
	RANDOMESTSETS="${FUNCTIONALOUTPUTDIR}/random_estimated_CNEs"
	bedtools shuffle -i "$ALLESTSETS" -g "$HGCHRINFOBEDTOOLS" -chrom -seed 66 > "$RANDOMESTSETS"
	head -n3 "$RANDOMESTSETS"
	wc -l "$RANDOMESTSETS"

	## OBTAIN DNASE AND H3K27AC SIGNAL FOR CNES AND RANDOM REGIONS

	for epigenetic_mark in "${ENCODESETS[@]}"; do
		output_dir="${FUNCTIONALOUTPUTDIR}/${epigenetic_mark}"
		mkdir -p "$output_dir"
		var_with_filename="ENCODESCREEN${epigenetic_mark}LIST"
		ENCODE_list="${!var_with_filename}"
		[[ -f "$ENCODE_list" ]] || error_exit "$LINENO" "Cannot find $ENCODE_list!"
		echo "# $ENCODE_list: $epigenetic_mark"

		while IFS= read -r line; 
		do
			name=$(echo $line | awk '{print $1}')
			echo "# $name"
			bigwig_file="$DEMOFILEDIR/"$(echo $line | awk '{print $3}')
			[[ -f "$bigwig_file" ]] || error_exit "$LINENO" "Cannot find $bigwig_file!"
			
			output="$output_dir/all_estimated_CNEs_${name}_mean"
			command="python ${SCRIPTS}/getbigwigsignalforbed_v1.py 8 $bigwig_file $ALLESTSETS > $output"		
			echo "# $command"
			eval "$command"

			output="$output_dir/random_estimated_CNEs_${name}_mean"
			command="python ${SCRIPTS}/getbigwigsignalforbed_v1.py 8 $bigwig_file $RANDOMESTSETS > $output"
			echo "# $command"
			eval "$command"

		done < <(grep -P -v "^#" "$ENCODE_list")
	done

fi


######################################################################################
# OVERLAPS WITH ENCODE SCREEN ELEMENTS
######################################################################################

# Overlap with cell-type agnostic SCREEN elements

if [[ "$OVERLAPS_WITH_CELLAGNOSTIC_SCREEN" = true ]]; then
	echo "# Overlap with cell-type agnostic SCREEN elements"

	output_dir="${FUNCTIONALOUTPUTDIR}/CellTypeAgnosticSCREEN"
	mkdir -p "$output_dir"

	input="$SUMSCREENLIST"
	{
		read 
		while IFS= read -r line
		do		
			for clade in "${CladesOfInterest[@]}"; do
				echo "# $clade"
				cne_file="${MAFOUTPUTDIR}/${clade}_estimated_CNEs_small"
				[[ -f "$cne_file" ]] || error_exit "$LINENO" "Cannot find $cne_file!"
				screen_file="$DEMOFILEDIR/"$(echo "$line" | awk '{print $2}')
				[[ -f "$screen_file" ]] || error_exit "$LINENO" "Cannot find $screen_file!"
				output=$(basename "$screen_file")
				output="${output_dir}/${clade}_estimated_CNEs_AND_${output%.bed}_SumFisher"
				cat "$screen_file" | sort -k1,1V -k2,2n > "$tmp1"
				command="bedtools fisher -a $cne_file -b $tmp1 -g $HGCHRINFOBEDTOOLS > $output"
				echo "# Running $command..."
				eval "$command"
				head "$output"
				wc -l "$output"
				rm "$tmp1"
			done
		done
	} < "$input"

	# Extract the relevant information from the output (p-value and ratio) and visualize
	output="${FUNCTIONALOUTPUTDIR}/fisher_sum_summary.txt"
	rm -f "$output"
	for file in "${output_dir}/"*_estimated_CNEs_AND*_SumFisher; do
		totalCNEs=$(cat $file | sed -n '1p' | cut -d' ' -f6)
		matchedCNEs=$(cat $file | sed -n '3p' | cut -d' ' -f5)
		SCREENelements=$(cat $file | sed -n '2p' | cut -d' ' -f6)
		number1=$(cat $file | tail -n 1 | awk '{print $3}')
		number2=$(cat $file | tail -n 1 | awk '{print $4}')
		echo -e "${file##*/}""\t"$matchedCNEs"\t"$totalCNEs"\t"$SCREENelements"\t"$number1"\t"$number2 
	done > "$output"
	head "$output"
	wc -l "$output"

	# Delete intermediate files
	rm -r "$output_dir"

	# Visualize
	command="R --slave --no-save --no-restore --no-environ --silent --args $output $FUNCTIONALOUTPUTDIR < ${SCRIPTS}/plot_sum_fisher_v2.R"
	echo "# $command"
	eval "$command"

fi

## Overlap with cell-type specific SCREEN elements

if [[ "$OVERLAPS_WITH_CELLSPEC_SCREEN" = true ]]; then
	echo "# Overlap with cell-type specific SCREEN elements"

	output_dir="${FUNCTIONALOUTPUTDIR}/CellTypeSpecificSCREEN"
	mkdir -p "$output_dir"

	input="$SINGLESCREENLIST"
	{
		read 
		while IFS= read -r line
		do		
			screen_file="$DEMOFILEDIR/"$(echo "$line" | awk '{print $4}')
			screen_tissue=$(echo "$line" | awk '{print $1}')
			[[ -f "$screen_file" ]] || error_exit "$LINENO" "Cannot find $screen_file!"
			cat "$screen_file" | sort -k1,1V -k2,2n > "$tmp1"
			
			screen_file=$(basename "$screen_file")

			for clade in "${CladesOfInterest[@]}";
			do
				cne_file="${MAFOUTPUTDIR}/${clade}_estimated_CNEs_small"
				[[ -f "$cne_file" ]] || error_exit "$LINENO" "Cannot find $cne_file!"
		
				for screen_set in $(cat "$tmp1" | cut -d$'\t' -f10 | sort | uniq);
				do
					echo "# $clade $screen_tissue $screen_set"

					output="${output_dir}/${clade}_estimated_CNEs_AND_${screen_tissue}_${screen_set}_SingleFisher"
					cat "$tmp1" | grep -P "\t$screen_set\t" > "$tmp2"
					wc -l "$tmp2"

					command="bedtools fisher -a $cne_file -b $tmp2 -g $HGCHRINFOBEDTOOLS > $output"
					echo "# Running $command..."
					eval "$command"
					head "$output"
					wc -l "$output"
					rm "$tmp2"
				done
			done
			rm "$tmp1"
		done
	} < "$input"

	# Extract the relevant information from the output (p-value and ratio) and visualize
	output="${FUNCTIONALOUTPUTDIR}/fisher_single_summary.txt"
	rm -f "$output"
	for file in "${output_dir}/"*_estimated_CNEs_AND*_SingleFisher;
	do
		totalCNEs=$(cat $file | sed -n '1p' | cut -d' ' -f6)
		matchedCNEs=$(cat $file | sed -n '3p' | cut -d' ' -f5)
		SCREENelements=$(cat $file | sed -n '2p' | cut -d' ' -f6)
		number1=$(cat $file | tail -n 1 | awk '{print $3}')
		number2=$(cat $file | tail -n 1 | awk '{print $4}')
		echo -e "${file##*/}""\t"$matchedCNEs"\t"$totalCNEs"\t"$SCREENelements"\t"$number1"\t"$number2 
	done > "$output"
	head "$output"
	wc -l "$output"

	# Delete intermediate files
	rm -r "$output_dir"

	# Visualize
	command="R --slave --no-save --no-restore --no-environ --silent --args $output $FUNCTIONALOUTPUTDIR < ${SCRIPTS}/plot_single_fisher_v2.R"
	echo "# $command"
	eval "$command"

fi

######################################################################################
# DNASE AND H3K27AC ANALYSIS
######################################################################################

# Perform DNase analysis

if [[ "$DNASE_ANALYSIS" = true ]]; then
	echo "# Performing DNase analysis..."

	epigenetic_mark="DNase"
	output_dir="${FUNCTIONALOUTPUTDIR}/${epigenetic_mark}"
	mkdir -p "$output_dir"

	grep -v '^#' "$ENCODESCREENDNaseLIST" | while IFS= read -r line; do
		tissue=$(echo $line | awk '{print $1}')
		echo "# $tissue"
		
		# Compute fold-changes relative to control random regions
		cne_mean="${output_dir}/all_estimated_CNEs_${tissue}_mean"
		random_mean="${output_dir}/random_estimated_CNEs_${tissue}_mean"
		command="R --slave --no-save --no-restore --no-environ --silent --args $cne_mean $random_mean $output_dir < ${SCRIPTS}/DNase_mean_to_fc_v2.R"
		echo "# $command"
		eval "$command"
		cne_fold_change="${output_dir}/all_estimated_CNEs_${tissue}_mean_corr"
		[[ -f "$cne_fold_change" ]] || error_exit "$LINENO" "Cannot find $cne_fold_change!"

		command="R --slave --no-save --no-restore --no-environ --silent --args $random_mean $random_mean $output_dir < ${SCRIPTS}/DNase_mean_to_fc_v2.R"
                echo "# $command"
                eval "$command"
		random_fold_change="${output_dir}/random_estimated_CNEs_${tissue}_mean_corr"
                [[ -f "$random_fold_change" ]] || error_exit "$LINENO" "Cannot find $random_fold_change!"

		# Assign categories based on the fold-changes
		command="R --slave --no-save --no-restore --no-environ --silent --args $cne_fold_change $output_dir < ${SCRIPTS}/add_category_fc_v2.R"
		echo "# $command"
		eval "$command"
	done

	# Visualize signal categories
	clade="all_estimated_CNEs_"
	targetstring="_mean_corr"
	command="R --slave --no-save --no-restore --no-environ --silent --args $clade $targetstring $SUM $output_dir < ${SCRIPTS}/EpigeneticCategory_analysis_v2.R"
	echo "# $command"
	eval "$command"

	# Enrichment calculation
	string2="random_estimated_CNEs_"
	command="R --slave --no-save --no-restore --no-environ --silent --args $clade $string2 $targetstring $SUM $output_dir < ${SCRIPTS}/Enrichment_v2.R"
	echo "# $command"
	eval "$command"

fi

# Perform H3k27ac analysis

if [[ "$H3K27AC_ANALYSIS" = true ]]; then
	echo "# Performing H3k27ac analysis..."
	
	epigenetic_mark="H3K27ac"
	output_dir="${FUNCTIONALOUTPUTDIR}/${epigenetic_mark}"
	mkdir -p "$output_dir"

	grep -v '^#' "$ENCODESCREENDNaseLIST" | while IFS= read -r line; do
		tissue=$(echo $line | awk '{print $1}')
		echo "# $tissue"

		# Compute fold-changes relative to control random regions
		cne_mean="${output_dir}/all_estimated_CNEs_${tissue}_mean"
		random_mean="${output_dir}/random_estimated_CNEs_${tissue}_mean"

		# Transform the means so that negative values become zeros
		command="R --slave --no-save --no-restore --no-environ --silent --args $cne_mean $output_dir < ${SCRIPTS}/correct_means_v2.R"
		echo "# $command"
		eval "$command"
		cne_fold_change="${output_dir}/all_estimated_CNEs_${tissue}_mean_corr"
		[[ -f "$cne_fold_change" ]] || error_exit "$LINENO" "Cannot find $cne_fold_change!"
		
		command="R --slave --no-save --no-restore --no-environ --silent --args $random_mean $output_dir < ${SCRIPTS}/correct_means_v2.R"
                echo "# $command"
                eval "$command"
		random_fold_change="${output_dir}/random_estimated_CNEs_${tissue}_mean_corr"
                [[ -f "$random_fold_change" ]] || error_exit "$LINENO" "Cannot find $random_fold_change!"

		# Assign categories based on the fold-changes
		command="R --slave --no-save --no-restore --no-environ --silent --args $cne_fold_change $output_dir < ${SCRIPTS}/add_category_fc_v2.R"
		echo "# $command"
		eval "$command"
	done

	# Visualize signal categories
	clade="all_estimated_CNEs_"
	targetstring="_mean_corr"
	command="R --slave --no-save --no-restore --no-environ --silent --args $clade $targetstring $SUM $output_dir < ${SCRIPTS}/EpigeneticCategory_analysis_v2.R"
	echo "# $command"
	eval "$command"

	# Enrichment calculation
	string2="random_estimated_CNEs_"
	targetstring="_mean_corr"
	command="R --slave --no-save --no-restore --no-environ --silent --args $clade $string2 $targetstring $SUM $output_dir < ${SCRIPTS}/Enrichment_v2.R"
	echo "# $command"
	eval "$command"
fi

# Summarize activity
ACTIVITY_ANALYSIS=true
if [[ "$ACTIVITY_ANALYSIS" = true ]]; then
	echo "# Consensus of DNase and H3K27ac categories"
	
	output_dir="$CONSDIR"
	mkdir -p "$output_dir"

	DNaseCAT="${FUNCTIONALOUTPUTDIR}/DNase/"$(basename "$ALLESTSETS")"__Categories.tsv"
	[[ -f "$DNaseCAT" ]] || error_exit "$LINENO" "Cannot find $DNaseCAT!"
	
	H3K27acCAT="${FUNCTIONALOUTPUTDIR}/H3K27ac/"$(basename "$ALLESTSETS")"__Categories.tsv"
	[[ -f "$H3K27acCAT" ]] || error_exit "$LINENO" "Cannot find $H3K27acCAT!"

	# Extract categories of interest
	command="R --slave --no-save --no-restore --no-environ --silent --args $DNaseCAT $H3K27acCAT $output_dir < ${SCRIPTS}/combine_categories_v2.R"
	echo "# $command"
	eval "$command"

	consensus="$CNE_ACTIVITY"
	[[ -f "$consensus" ]] || error_exit "$LINENO" "Cannot find $consensus!"
	command="R --slave --no-save --no-restore --no-environ --silent --args $consensus ${SUM} $output_dir < ${SCRIPTS}/Consensus_analysis_v3.R"
	echo "# $command"
	eval "$command"
fi

######################################################################################
# OVERLAP WITH PROMOTERS
######################################################################################

if [[ "$OVERLAPS_WITH_PROMOTERS" = true ]]; then
	echo "# Computing overlap with promoters..."

	for clade in "${CladesOfInterest[@]}"; do
		cne_file="${MAFOUTPUTDIR}/${clade}_estimated_CNEs_small"
		[[ -f "$cne_file" ]] || error_exit "$LINENO" "Cannot find $cne_file!"

		command="bedtools intersect -a $cne_file -b $TSSRELEVANT -wao > $tmp1" 
		echo "# $command"
		eval "$command"
		head -n3 "$tmp1"
		wc -l "$tmp1"
		
		command="R --slave --no-save --no-restore --no-environ --silent --args ${clade} ${tmp1} ${FUNCTIONALOUTPUTDIR} < ${SCRIPTS}/get_CNE_w-wo_TSS_promoter_v2.R"
		echo "# $command"
		eval "$command"
		rm "$tmp1"
	done
fi

echo "Done!"
