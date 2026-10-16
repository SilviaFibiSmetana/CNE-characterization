#!/bin/bash

# Input Data
[[ -f "$GENOMEGFF" ]] || error_exit "$LINENO" "Cannot find $GENOMEGFF!"

[[ -f "$PHASTCONSELEMENTSSELECTEDLENGTHFILE" ]] || error_exit "$LINENO" "Cannot find $PHASTCONSELEMENTSSELECTEDLENGTHFILE!"

# File prefixes

GENOMEGFFNAME="${GENOMEGFF##*/}"

# Output

[[ -d "$ANNOTATIONOUTPUTDIR" ]] || error_exit "$LINENO" "Cannot find $ANNOTATIONOUTPUTDIR!"

ELEMENTSCODING="${ANNOTATIONOUTPUTDIR}/coding"
ELEMENTSNONCODING="${ANNOTATIONOUTPUTDIR}/noncoding"

ELEMENTSNONCODINGINTROINTER="${ANNOTATIONOUTPUTDIR}/noncoding_introinter.txt"
ELEMENTSNONCODINGINTRO="${ANNOTATIONOUTPUTDIR}/noncoding_intro.txt"
ELEMENTSNONCODINGINTER="${ANNOTATIONOUTPUTDIR}/noncoding_inter.txt"

# Temporary files

tmp1=$(mktemp -u -p . --suffix=".bed")
tmp2=$(mktemp -u -p .)
tmp3=$(mktemp -u -p .)
tmp4=$(mktemp -u -p .)

# Get only the lines of the coding sequences cds and extract the important colums and parts of the columns(ens-id) and sort them

zcat "$GENOMEGFF" | awk '{if($3 == "CDS"){printf "%s\t%s\t%s\t%s\n", $1,$4,$5,$9}}' | perl -ne '{if(/^([^\s]+\s\d+\s\d+)\s[^\s]+gene_id=(ENS.[\d.]+)/){print "$1\t$2\n";}}' | sort -k1,1 -k2,2n > "$tmp1"

# Remove the duplicates of the sorted list of cds and merge the cds which overlap or are direct neighbors

GENOMEGFFCODING="${ANNOTATIONOUTPUTDIR}/${GENOMEGFFNAME%.gff3.gz}_CDS.bed"

bedtools merge -c 4 -o distinct -delim ";" -i "$tmp1" | uniq > "$GENOMEGFFCODING"
head -n3 "$GENOMEGFFCODING"
wc -l "$GENOMEGFFCODING"

rm $tmp1

echo

# Get only the rna regions of the genome and convert to bed format

GENOMEGFFRNA="${ANNOTATIONOUTPUTDIR}/${GENOMEGFFNAME%.gff3.gz}_RNA.bed"
		
zcat "$GENOMEGFF" | awk '{printf "%s\t%s\t%s\t%s\t%s\n", $1,$3,$4,$5,$9}' | perl -ne '{if(/^([^\s]+\s[^\s]+\s\d+\s\d+)\s[^\s]*gene_type=([^\s;]+)/){print "$1\t$2\n";}}' | grep "RNA" | awk '{if($2=="exon") print $1"\t"$3"\t"$4"\t"$5}' | sort -k1,1 -k2,2n > "$tmp1";
head -n3 "$tmp1"
wc -l "$tmp1"

# Remove the duplicates of the sorted list of cds and merge the cds which overlap or are direct neighbors

bedtools merge -c 4 -o distinct -delim ";" -i "$tmp1" | uniq > "$GENOMEGFFRNA"
rm "$tmp1"
head -n3 "$GENOMEGFFRNA"
wc -l "$GENOMEGFFRNA"

echo

# Intersect phastcons with cds regions and decide whether the phastcon is coding or noncoding

INTERSECTIONRAWCDS="${ANNOTATIONOUTPUTDIR}/intersection_rawCDS.txt"
wc -l "$PHASTCONSELEMENTSSELECTEDLENGTHFILE"
bedtools intersect -wao -a "$PHASTCONSELEMENTSSELECTEDLENGTHFILE" -b "$GENOMEGFFCODING" > "$INTERSECTIONRAWCDS"
head -n5 "$INTERSECTIONRAWCDS"
wc -l "$INTERSECTIONRAWCDS"

cut -f1-3 "$INTERSECTIONRAWCDS" | sort -k1,1 -k2,2n | uniq | wc -l
head -n5 "$INTERSECTIONRAWCDS"
wc -l "$INTERSECTIONRAWCDS"

echo

# Intersect phastcons with RNA regions

INTERSECTIONRAWRNA="${ANNOTATIONOUTPUTDIR}/intersection_rawRNA.txt"
head "$PHASTCONSELEMENTSSELECTEDLENGTHFILE"
wc -l "$PHASTCONSELEMENTSSELECTEDLENGTHFILE"
bedtools intersect -wao -a "$PHASTCONSELEMENTSSELECTEDLENGTHFILE" -b "$GENOMEGFFRNA" > "$INTERSECTIONRAWRNA"
head -n5 "$INTERSECTIONRAWRNA"
wc -l "$INTERSECTIONRAWRNA"

cut -f1-3 "$INTERSECTIONRAWRNA" | sort -k1,1 -k2,2n | uniq | wc -l
head -n5 "$INTERSECTIONRAWRNA"
wc -l "$INTERSECTIONRAWRNA"

echo

# Find phastcons overlapping with multiple cds, therefore first combine the content of the first 4 columns

INTERSECTIONRAWCDS1="${ANNOTATIONOUTPUTDIR}/intersection_rawCDS1.txt"
INTERSECTIONANNCDS="${INTERSECTIONRAWCDS1}_annotated"
ELEMENTSSTILLPHASTCON="${ANNOTATIONOUTPUTDIR}/stillphastcon"
cat "$INTERSECTIONRAWCDS" | awk '{print $1"_"$2"_"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' > "$tmp1"
head -n5 "$tmp1"
wc -l "$tmp1"

echo

cat "$tmp1" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > "$INTERSECTIONRAWCDS1"
rm "$tmp1"
head "$INTERSECTIONRAWCDS1"
wc -l "$INTERSECTIONRAWCDS1"

echo

R --slave --no-save --no-restore --no-environ --silent --args "$INTERSECTIONRAWCDS1" "$ANNOTATIONOUTPUTDIR" < "${SCRIPTS}/phastovmult_basic_annotation.R"
head -n5 "$INTERSECTIONANNCDS"
wc -l "$INTERSECTIONANNCDS"
cat "$INTERSECTIONANNCDS" | grep -v "noncoding" > "$ELEMENTSCODING"
cat "$INTERSECTIONANNCDS" | grep "noncoding" > "$ELEMENTSSTILLPHASTCON"
wc -l "$ELEMENTSCODING"
wc -l "$ELEMENTSSTILLPHASTCON"

echo

# For the noncoding elements, find phastcons overlapping with multiple rna, therefore first combine the content of the first 4 columns
INTERSECTIONRAWRNA1="${ANNOTATIONOUTPUTDIR}/intersection_rawRNA1.txt"
INTERSECTIONANNRNA="${INTERSECTIONRAWRNA1}_annotated"
ELEMENTSRNA="${ANNOTATIONOUTPUTDIR}/RNA"

cat "$ELEMENTSSTILLPHASTCON" | awk '{print $1"_"$2"_"$3}' | sort -k1,1 > "$tmp1"
head "$tmp1"
wc -l "$tmp1";

cat "$INTERSECTIONRAWRNA" | awk '{print $1"_"$2"_"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' > "$tmp2"
head -n5 "$tmp2"
wc -l "$tmp2"

echo

cat "$tmp2" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > "$tmp3"
join -1 1 -2 1 -t $'\t' <(sort -k1,1 "$tmp3") <(sort -k1,1 "$tmp1") > "$tmp4"
cat "$tmp4" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > "$INTERSECTIONRAWRNA1"
rm "$tmp1" "$tmp2" "$tmp3" "$tmp4"
head "$INTERSECTIONRAWRNA1"
wc -l "$INTERSECTIONRAWRNA1"

echo

R --slave --no-save --no-restore --no-environ --silent --args "$INTERSECTIONRAWRNA1" "$ANNOTATIONOUTPUTDIR" < "${SCRIPTS}/phastovmult_basic_annotation.R";
head -n5 "$INTERSECTIONANNRNA"
wc -l "$INTERSECTIONANNRNA"

cat "$INTERSECTIONANNRNA" | grep -v "noncoding" > "$ELEMENTSRNA"
wc -l "$ELEMENTSRNA"

cat "$INTERSECTIONANNRNA" | grep "noncoding" > "$ELEMENTSNONCODING"
wc -l "$ELEMENTSNONCODING"

echo

# Define intronic and intergenic cnes by overlapping them with the gene lines of the gff file (as we made already before sure they are not coding)

GENOMEGFFGENE="${ANNOTATIONOUTPUTDIR}/${GENOMEGFFNAME%.gff3.gz}_gene.bed"

zcat "$GENOMEGFF" | awk '{if($3=="gene") print $1"\t"$4"\t"$5}' > "${GENOMEGFFGENE}"
wc -l "$GENOMEGFFGENE"

cat "$ELEMENTSNONCODING" | sort -k1,1V -k2,2n > "$tmp1"
wc -l "$tmp1"

bedtools intersect -a "$tmp1" -b "${GENOMEGFFGENE}" -wao > "${ANNOTATIONOUTPUTDIR}/CNEs_gff_gene_intersected.txt"
rm "$tmp1"
head "${ANNOTATIONOUTPUTDIR}/CNEs_gff_gene_intersected.txt"
wc -l "${ANNOTATIONOUTPUTDIR}/CNEs_gff_gene_intersected.txt"

echo

R --slave --no-save --no-restore --no-environ --silent --args "${ANNOTATIONOUTPUTDIR}/CNEs_gff_gene_intersected.txt" "$ANNOTATIONOUTPUTDIR" < "${SCRIPTS}/phastovmult_cne_annotation.R"

echo "# Annotated conserved elements are stored in"

head -n3 "$ELEMENTSNONCODINGINTROINTER"
wc -l "$ELEMENTSNONCODINGINTROINTER"
head -n3 "$ELEMENTSNONCODINGINTRO"
wc -l "$ELEMENTSNONCODINGINTRO"
head -n3 "$ELEMENTSNONCODINGINTER"
wc -l "$ELEMENTSNONCODINGINTER"

echo "Done!"
