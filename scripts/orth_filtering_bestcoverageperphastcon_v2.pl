#! /usr/bin/perl

use strict;
use warnings;

##usage: perl orth_filtering_bestcoverageperphastcon.pl file_orthologs_linedup outputfile

##get files from commandline
unless($#ARGV == 1){
	die "usage: perl $0 file_orthologs_linedup outputfile !\n";
}

my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];


##define variables 

my $phastchr1;
my $phaststart1;
my $phastend1;
my $hgchr1;
my $hgstrand1;
my $hgdist1;
my $orthgenome1;
my $orthchr1;
my $orthstart1;
my $orthend1;
my $orthseq1;
my $orthstrand1;
my $orthlength1;
my $coverage1;

my $phastchr2;
my $phaststart2;
my $phastend2;
my $hgchr2;
my $hgstrand2;
my $hgdist2;
my $orthgenome2;
my $orthchr2;
my $orthstart2;
my $orthend2;
my $orthseq2;
my $orthstrand2; 
my $orthlength2;
my $coverage2;

my $eofmarker;

##read and process file_orthologs_merged_overlapcorr
open(READ, $inputfile) || die "Cannot open $inputfile!\n";
open(OUT, '>', $outputfile) or die $!;

my $first_line_done = 0;
while(<READ>){
	chomp;
	##process the first line separately to have this line to compare it with the second
	if(!$first_line_done) {
		if( /^([^\s]+)\s(\d+)\s(\d+)\s([^\s]+)\s([+-])\s([^\s]+)\s([^\s]+)\s([^\s]+)\s(\d+)\s(\d+)\s([ACTGNXactgnx-]+)\s([+-])\s(\d+)\s([\d.]+)/) {
			#print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$12\t$13\t$14\n";
			
			$phastchr1 = $1;
			$phaststart1 = $2;
			$phastend1 = $3;
			$hgchr1 = $4;
			$hgstrand1 = $5;
			$hgdist1 = $6;
			$orthgenome1 = $7;
			$orthchr1 = $8;
			$orthstart1 = $9;
			$orthend1 = $10;
			$orthseq1 = $11;
			$orthstrand1 = $12;
			$orthlength1 = $13;
			$coverage1 = $14;	
			
			$first_line_done++;
		}else{
			die "$_: Wrong format of first line!\n";
		}
	}else{
	##process all the other lines
		if( /^([^\s]+)\s(\d+)\s(\d+)\s([^\s]+)\s([+-])\s([^\s]+)\s([^\s]+)\s([^\s]+)\s(\d+)\s(\d+)\s([ACTGNXactgnx-]+)\s([+-])\s(\d+)\s([\d.]+)/) {
			#print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$12\t$13\t$14\n";
			
			$phastchr2 = $1;
			$phaststart2 = $2;
			$phastend2 = $3;
			$hgchr2 = $4;
			$hgstrand2 = $5;
			$hgdist2 = $6;
			$orthgenome2 = $7;
			$orthchr2 = $8;
			$orthstart2 = $9;
			$orthend2 = $10;
			$orthseq2 = $11;
			$orthstrand2 = $12;
			$orthlength2 = $13;
			$coverage2 = $14;
			
			
			#filter for the ortholog with the highest coverage
			if( ($phastchr2 eq $phastchr1) && ($phaststart2 eq $phaststart1) && ($phastend2 eq $phastend1) && ($hgchr2 eq $hgchr1) && ($hgstrand2 eq $hgstrand1) && ($orthgenome2 eq $orthgenome1) ){
			
				if ($coverage1 gt $coverage2) {
					#I keep line 1
					$eofmarker = 0;
				}else{
					#I keep line 2 and move it into line 1
					$phastchr1 = $phastchr2;
					$phaststart1 = $phaststart2;
					$phastend1 = $phastend2;
					$hgchr1 = $hgchr2;
					$hgstrand1 = $hgstrand2;
					$hgdist1 = $hgdist2;
					$orthgenome1 = $orthgenome2;
					$orthchr1 = $orthchr2;
					$orthstart1 = $orthstart2;
					$orthend1 = $orthend2;
					$orthseq1 = $orthseq2;
					$orthstrand1 = $orthstrand2;
					$orthlength1 = $orthlength2;
					$coverage1 = $coverage2;
					
					$eofmarker = 1;
				}
			
			}else{
				print OUT "$phastchr1\t$phaststart1\t$phastend1\t$hgchr1\t$hgstrand1\t$hgdist1\t$orthgenome1\t$orthchr1\t$orthstart1\t$orthend1\t$orthseq1\t$orthstrand1\t$orthlength1\t$coverage1\n";

			#zeile mit höchster coverage printen
			#alle 2 zu 1 um sie mit nächster zeile auszugeben
				$phastchr1 = $phastchr2;
				$phaststart1 = $phaststart2;
				$phastend1 = $phastend2;
				$hgchr1 = $hgchr2;
				$hgstrand1 = $hgstrand2;
				$hgdist1 = $hgdist2;
				$orthgenome1 = $orthgenome2;
				$orthchr1 = $orthchr2;
				$orthstart1 = $orthstart2;
				$orthend1 = $orthend2;
				$orthseq1 = $orthseq2;
				$orthstrand1 = $orthstrand2;
				$orthlength1 = $orthlength2;
				$coverage1 = $coverage2;
				
				$eofmarker = 0;
			}
			
			if(eof(READ) && ($eofmarker == 0)){
				print OUT "$phastchr1\t$phaststart1\t$phastend1\t$hgchr1\t$hgstrand1\t$hgdist1\t$orthgenome1\t$orthchr1\t$orthstart1\t$orthend1\t$orthseq1\t$orthstrand1\t$orthlength1\t$coverage1\n";
			}elsif(eof(READ) && ($eofmarker == 1)){
				print OUT "$phastchr2\t$phaststart2\t$phastend2\t$hgchr2\t$hgstrand2\t$hgdist2\t$orthgenome2\t$orthchr2\t$orthstart2\t$orthend2\t$orthseq2\t$orthstrand2\t$orthlength2\t$coverage2\n";
			
			}
			
		}else{
			die "$_: Wrong format of line!\n";
		}
	}
}

close(READ);
close(OUT);

exit 0;

