#!/bin/bash

set -Eeuo pipefail

######################################################
##
## This script downloads and preprocesses
## datasets from the ENCODE data repository.
##
#######################################################

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

######################################################################################
# DATASETS
######################################################################################

ENCODESCREENDNaseLIST="ENCODE_SCREEN_list_DNase.txt"
ENCODESCREENH3K27acLIST="ENCODE_SCREEN_list_H3K27ac.txt"
[[ -f "$ENCODESCREENDNaseLIST" ]] || error_exit "$LINENO" "Cannot find $ENCODESCREENDNaseLIST!"
[[ -f "$ENCODESCREENH3K27acLIST" ]] || error_exit "$LINENO" "Cannot find $ENCODESCREENH3K27acLIST!"

DNASE=true
H3K27AC=true

######################################################################################
# DOWNLOAD BAM FILES AND TRANSFORM TO BIGWIG
######################################################################################

# Selected ENCODE DNAse Datasets

if [[ "$DNASE" = true ]]; then
	echo "# Downloading ENCODE DNase Datasets in $ENCODESCREENDNaseLIST..."

	grep -v '^#' "$ENCODESCREENDNaseLIST" | while IFS= read -r line; do
		tissue=$(echo $line | awk '{print $1}')
		ENCODEId=$(echo $line | awk '{print $2}')
		bigWig_file=$(echo $line | awk '{print $3}')

		echo "# $tissue: $ENCODEId and $control"
		# bam file
		target=$(echo "$bigWig_file" | sed "s/.bigWig$//")".bam"
		output_dir=$(dirname "$target")
		[[ -d "$output_dir" ]] || mkdir -p "$output_dir"
		
		# download bam file
		command="wget -O $target https://www.encodeproject.org/files/${ENCODEId}/@@download/${ENCODEId}.bam"
		echo "# $command"
		eval "$command"
		# index bam file
		command="samtools index $target"
		echo "# $command"
		eval "$command"
		# compute coverage/convert to bigWig format
		command="bamCoverage --binSize 10 --extendReads 300 --normalizeUsing RPKM -p 8 --bam $target --outFileFormat bigwig --outFileName $bigWig_file"
		echo "# $command"
		eval "$command"
		# delete bam file
		rm "$target"
	done
fi

# Selected ENCODE H3K27ac Datasets

if [[ "$H3K27AC" = true ]]; then
	echo "# Downloading ENCODE H3k27ac Datasets in $ENCODESCREENH3K27acLIST..."

	grep -v '^#' "$ENCODESCREENH3K27acLIST" | while IFS= read -r line; do
		tissue=$(echo $line | awk '{print $1}')
		ENCODEIds=$(echo $line | awk '{print $2}')
		bigWig_file=$(echo $line | awk '{print $3}')
		control=$(echo $line | awk '{print $4}')
		
		echo "# $tissue: $ENCODEIds and $control"

		# download bam files
		output_dir=$(dirname "$bigWig_file")
		[[ -d "$output_dir" ]] || mkdir -p "$output_dir"
		tr ',' '\n' <<< "$ENCODEIds" | while read -r ENCODEId; do
			target="$output_dir/${ENCODEId}.bam"
			command="wget -O $target https://www.encodeproject.org/files/${ENCODEId}/@@download/${ENCODEId}.bam"
			echo "# $command"
			eval "$command"
		done
		target="$output_dir/${control}.bam"
		command="wget -O $target $output_dir https://www.encodeproject.org/files/${control}/@@download/${control}.bam"
		echo "# $command"
		eval "$command"
	
		# merge bam files
		target="$output_dir/${tissue}.bam"
		command="samtools merge -o $target $output_dir/{$ENCODEIds}.bam"
		echo "# $command"
		eval "$command"

		# index bam files
		target="$output_dir/${tissue}.bam"
		command="samtools index $target"
		echo "# $command"
		eval "$command"

		target="$output_dir/${control}.bam"
		command="samtools index $target"
		echo "# $command"
		eval "$command"

		# compute coverage/convert to bigWig format
		command="bamCompare --binSize 10 --extendReads 300 --scaleFactorsMethod SES -p 8 --skipNonCoveredRegions -b1 $output_dir/${tissue}.bam -b2 $output_dir/${control}.bam --outFileFormat bigwig --outFileName $bigWig_file"
		echo "# $command"
		eval "$command"
		
		# delete bam files
		command="rm $output_dir/{$ENCODEIds}.bam*"
		eval "$command"
		command="rm $output_dir/${control}.bam*"
		eval "$command"
	done
fi

echo "Done!"

