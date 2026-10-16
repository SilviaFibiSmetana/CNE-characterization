#! /usr/bin/perl

use strict;
use warnings;

##usage: perl merge_to_phastcon.pl file_orthextracted_strandcorrectedsorted outputfile allowed_distance_as_number

##get files from commandline
unless($#ARGV == 2){
	die "usage: perl $0 file_orthextracted_strandcorrectedsorted outputfile allowed_distance_as_number!\n";
}

my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];
my $alloweddist = $ARGV[2];

##define variables 

my $phastchr1;
my $phaststart1;
my $phastend1;
my $hgchr1;
my $hgstart1;
my $hgend1;
my $hgstrand1;
my $orthgenome1;
my $orthchr1;
my $orthstart1;
my $orthend1;
my $orthseq1;
my $orthstrand1;
my $score1;
my $length1;
my $orthlength1;

my $phastchr2;
my $phaststart2;
my $phastend2;
my $hgchr2;
my $hgstart2;
my $hgend2;
my $hgstrand2;
my $orthgenome2;
my $orthchr2;
my $orthstart2;
my $orthend2;
my $orthseq2;
my $orthstrand2;
my $score2;
my $length2;
my $orthlength2;

my $hglength1;
my $hglength2;
my $coverage1;
my $coverage2;
my $hgstart1mod;
my $hgdistcalc = 0;
my $hgdist = "noorth";

##read and process file_orthologs_merged_overlapcorr
open(READ, $inputfile) || die "Cannot open $inputfile!\n";
open(OUT, '>', $outputfile) or die $!;

my $first_line_done = 0;
while(<READ>){
	chomp;
	##process the first line separately to have this line to compare it with the second
	if(!$first_line_done) {
		if( /([^\s]+)\s+(\d+)\s+(\d+)\s+([^\s]+)\s+(\d+)\s+(\d+)\s+([+-])\s+([^\s]+)\s+([^\s]+)\s+(\d+)\s+(\d+)\s+([ACTGNactgn-]+)\s+([+-])\s+([^\s]+)\s+(\d+)\s+(\d+)/ ) {
			#print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$12\t$13\t$14\t$15\t$16\n";
			
			$phastchr1 = $1;
			$phaststart1 = $2;
			$phastend1 = $3;
			$hgchr1 = $4;
			$hgstart1 = $5;
			$hgend1 = $6;
			$hgstrand1 = $7;
			$orthgenome1 = $8;
			$orthchr1 = $9;
			$orthstart1 = $10;
			$orthend1 = $11;
			$orthseq1 = $12;
			$orthstrand1 = $13;
			$score1 = $14;
			$length1 = $15;
			$orthlength1 = $16;
			
			$hgstart1mod = $hgstart1;
			$hglength1 = ($hgend1 - $hgstart1 + 1);
			$coverage1 = $hglength1 / ($phastend1 - $phaststart1 + 1);
			
			$first_line_done++;
		}else{
			die "$_: Wrong format of first line!\n";
		}
	}else{
	##process all the other lines
		if( /([^\s]+)\s+(\d+)\s+(\d+)\s+([^\s]+)\s+(\d+)\s+(\d+)\s+([+-])\s+([^\s]+)\s+([^\s]+)\s+(\d+)\s+(\d+)\s+([ACTGNactgn-]+)\s+([+-])\s+([^\s]+)\s+([^\s]+)\s+(\d+)/ ) {
			#print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$12\t$13\t$14\t$15\t$16\n";

			$phastchr2 = $1;
			$phaststart2 = $2;
			$phastend2 = $3;
			$hgchr2 = $4;
			$hgstart2 = $5;
			$hgend2 = $6;
			$hgstrand2 = $7;
			$orthgenome2 = $8;
			$orthchr2 = $9;
			$orthstart2 = $10;
			$orthend2 = $11;
			$orthseq2 = $12;
			$orthstrand2 = $13;
			$score2 = $14;
			$length2 = $15;
			$orthlength2 = $16;
			
			$hglength2 = ($hgend2 - $hgstart2 + 1);
			$coverage2 = $hglength2 / ($phastend2 - $phaststart2 + 1);

			##if the conditions match to line up two orthologs closer than $alloweddist		
			if( ($phastchr2 eq $phastchr1) && ($phaststart2 eq $phaststart1) && ($phastend2 eq $phastend1) && ($hgchr2 eq $hgchr1) && ($hgstrand2 eq $hgstrand1) && ($orthgenome2 eq $orthgenome1) && ($orthchr2 eq $orthchr1) && ($orthstrand2 eq $orthstrand1) && ($orthstart2 <= $orthend1+$alloweddist) ){
				
				$hglength1 = $hglength1 + $hglength2;
				$coverage1 = $hglength1 / ($phastend1 - $phaststart1 + 1);
				
				$phastchr1 = $phastchr2;
				$phaststart1 = $phaststart2;
				$phastend1 = $phastend2;
				$hgchr1 = $hgchr2;
				#$hgstart1 = $hgstart1;
				#$hgend1 = $hgend2;
				$hgstrand1 = $hgstrand2;
				$orthgenome1 = $orthgenome2;
				$orthchr1 = $orthchr2;
				#$orthstart1 = $orthstart1;
				$orthend1 = $orthend2;
				$orthseq1 = $orthseq1.'x'.$orthseq2;
				$orthstrand1 = $orthstrand2;
				#$score1 = $score1;
				$length1 = $length1+$length2;
				$orthlength1 = $orthlength1+$orthlength2;
				
				$hgdistcalc = $hgstart2 - $hgend1;
				
				if( ($orthstrand1 eq "-")){
					$hgdistcalc = $hgstart1mod - $hgend2;
				}elsif( ($orthstrand1 eq "+")){
					$hgdistcalc = $hgstart2 - $hgend1;
				}
				
				if( ($hgdistcalc lt 0)){
					$hgdistcalc = -1;
				}
				
				$hgdist = $hgdist.'_'.$hgdistcalc;
				$hgstart1mod = $hgstart2;
				
			}else{
				print OUT "$phastchr1\t$phaststart1\t$phastend1\t$hgchr1\t$hgstrand1\t$hgdist\t$orthgenome1\t$orthchr1\t$orthstart1\t$orthend1\t$orthseq1\t$orthstrand1\t$orthlength1\t$coverage1\n";
				
				$hglength1 = $hglength2;
				$coverage1 = $coverage2;
				
				$phastchr1 = $phastchr2;
				$phaststart1 = $phaststart2;
				$phastend1 = $phastend2;
				$hgchr1 = $hgchr2;
				$hgstart1 = $hgstart2;
				$hgend1 = $hgend2;
				$hgstrand1 = $hgstrand2;
				$orthgenome1 = $orthgenome2;
				$orthchr1 = $orthchr2;
				$orthstart1 = $orthstart2;
				$orthend1 = $orthend2;
				$orthseq1 = $orthseq2;
				$orthstrand1 = $orthstrand2;
				$score1 = $score2;
				$length1 = $length2;
				$orthlength1 = $orthlength2;
				$hgdist = "noorth";
				$hgstart1mod = $hgstart2;
			}

		}else{	
			die "$_: Wrong format!\n";
		}
	}
}

close(READ);
close(OUT);

