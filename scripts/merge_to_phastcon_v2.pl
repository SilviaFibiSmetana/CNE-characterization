#! /usr/bin/perl

use strict;
use warnings;

##usage: perl merge_to_phastcon.pl file_orthextracted_strandcorrectedsorted outputfile

##get files from commandline
unless($#ARGV == 1){
	die "usage: perl $0 file_orthextracted_strandcorrectedsorted outputfile!\n";
}

my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];

##define variables 
#my @line;

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
my $orthseq1;
my $orthstrand1;
my $score1;
my $length1;

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
my $orthseq2;
my $orthstrand2;
my $score2;
my $length2;

my $pureseq1;
my $pureseqlength1;
my $pureseq2;
my $pureseqlength2;

my $orthend1;
my $orthend2;

my $endoforth;

my $tmp;
my $orthlength;

##read and process file_orthextracted_strandcorrectedsorted
open(READ, $inputfile) || die "Cannot open $inputfile!\n";
open(OUT, '>', $outputfile) or die $!;

my $first_line_done = 0;
while(<READ>){
	chomp;
#process the first line separately to have this line to compare it with the second
	if(!$first_line_done) {
	#chr1	183128	183222	chr1	183128	183128	+	ailMel1	GL192379.1	2859287	C	+	183382	-1: Wrong format!

		if( /([^\s]+)\s+(\d+)\s+(\d+)\s+([^\s]+)\s+(\d+)\s+(\d+)\s+([+-])\s+([^\s]+)\s+([^\s]+)\s+(\d+)\s+([ACTGNactgn-]+)\s+([+-])\s+([^\s]+)\s+([\d+-]+)/ ) {
			#print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$12\t$13\t$14\n";
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
			$orthseq1 = $11;
			$orthstrand1 = $12;
			$score1 = $13;
			$length1 = $14;
			
			$pureseq1 = $orthseq1;
			$pureseq1 =~ s/-//g;
			$pureseqlength1 = length($pureseq1);
			$orthseq1 = $pureseq1;				
			$first_line_done++;
		}else{
			print "Wrong format of first line!\n";
		}
	}else{
#process all the other lines
		if( /([^\s]+)\s+(\d+)\s+(\d+)\s+([^\s]+)\s+(\d+)\s+(\d+)\s+([+-])\s+([^\s]+)\s+([^\s]+)\s+(\d+)\s+([ACTGNactgn-]+)\s+([+-])\s+([^\s]+)\s+([\d+-]+)/ ) {
			#print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$12\t$13\t$14\n";
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
			$orthseq2 = $11;
			$orthstrand2 = $12;
			$score2 = $13;
			$length2 = $14;
			
			$pureseq2 = $orthseq2;
			$pureseq2 =~ s/-//g;
			$pureseqlength2 = length($pureseq2);
			$orthseq2 = $pureseq2;

			if( ($orthstrand1 eq "-")){
				$orthend1 = $orthstart1 - $pureseqlength1 + 1;
			}elsif( ($orthstrand1 eq "+")){
				$orthend1 = $orthstart1 + $pureseqlength1 - 1;
			}
			$endoforth = 0;
#if the conditions match to combine two lines into one because they are the same ortholog		
			if( ($phastchr2 eq $phastchr1) && ($phaststart2 eq $phaststart1) && ($phastend2 eq $phastend1) && ($hgchr2 eq $hgchr1) && ($hgstart2 == $hgend1+1) && ($hgstrand2 eq $hgstrand1) && ($orthgenome2 eq $orthgenome1) && ($orthchr2 eq $orthchr1) && (($orthstart2 == $orthstart1+$pureseqlength1) || ($orthstart2 == $orthstart1-$pureseqlength1)) ){
				$phastchr1 = $phastchr2;
				$phaststart1 = $phaststart2;
				$phastend1 = $phastend2;
				$hgchr1 = $hgchr2;
				$hgend1 = $hgend2;
				$hgstrand1 = $hgstrand2;
				$orthgenome1 = $orthgenome2;
				$orthchr1 = $orthchr2;
				
				$orthseq1 = $orthseq1 . $orthseq2;
				$orthstrand1 = $orthstrand2;
				$score1 = $score2;
				$length1 = $length1+$length2;
				
				$pureseqlength1 = $pureseqlength1+$pureseqlength2;
				if( ($orthstrand1 eq "-")){ #$orthstart1 > $orthstart2 ){
					$orthend1 = $orthstart1 - $pureseqlength1 + 1;
				}elsif( ($orthstrand1 eq "+")){ #$orthstart1 < $orthstart2 ){
					$orthend1 = $orthstart1 + $pureseqlength1 - 1;
				}
				
#if the conditions match and we reach the end of the file, print the current ortholog			
				
			}else{
				#what to do when one ortholog is finished
				$endoforth = 1;
				if( $orthstrand1 eq "-" ){
					$tmp = $orthend1;
					$orthend1 = $orthstart1;
					$orthstart1 = $tmp;
				}
				$orthlength = $orthend1 - $orthstart1 + 1;

				if ( $orthstart1 > $orthend1 || $hgstart1 > $hgend1 ){
					die "Something went wrong! $!\n";
				}
				print OUT "$phastchr1\t$phaststart1\t$phastend1\t$hgchr1\t$hgstart1\t$hgend1\t$hgstrand1\t$orthgenome1\t$orthchr1\t$orthstart1\t$orthend1\t$orthseq1\t$orthstrand1\t$score1\t$length1\t$orthlength\n";

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
				$orthseq1 = $orthseq2;
				$orthstrand1 = $orthstrand2;
				$score1 = $score2;
				$length1 = $length2;
					
				$pureseqlength1 = $pureseqlength2;
			}
			if(eof(READ) && ($endoforth == 0)){
					if($orthstrand1 eq "-"){
						$tmp = $orthend1;
						$orthend1 = $orthstart1;
						$orthstart1 = $tmp;
					}
					$orthlength = $orthend1 - $orthstart1 + 1;
					
					print OUT "$phastchr1\t$phaststart1\t$phastend1\t$hgchr1\t$hgstart1\t$hgend1\t$hgstrand1\t$orthgenome1\t$orthchr1\t$orthstart1\t$orthend1\t$orthseq1\t$orthstrand1\t$score1\t$length1\t$orthlength\n";
				}
		}#elsif( /([!^\s]+)\s+(!\d+)\s+(!\d+)\s+(![^\s]+)\s+(!\d+)\s+(!\d+)\s+(![+-])\s+(![^\s]+)\s+(![^\s]+)\s+(!\d+)\s+(![ACTGNactgn-]+)\s+(![+-])\s+(![^ ]+)\s+(!\d+)/ ){
		else{	
			#print "Wrong Format!\n";
			die "$_: Wrong format!\n";
		}
	}
}

close(READ);
close(OUT);

exit 0;



	
