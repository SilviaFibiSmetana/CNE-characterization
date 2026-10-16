#! /usr/bin/perl

use strict;
use warnings;

##usage: perl fit_orthologs.pl file_orthologs_merged outputfile

##get files from commandline
unless($#ARGV == 1){
	die "usage: perl $0 file_orthologs_merged outputfile!\n";
}

my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];

##define variables 

my $phastchr;
my $phaststart;
my $phastend;
my $hgchr;
my $hgstart;
my $hgend;
my $hgstrand;
my $orthgenome;
my $orthchr;
my $orthstart;
my $orthend;
my $orthseq;
my $orthstrand;
my $score;
my $length;
my $reallength;

my $difference;


##read and process file_orthextracted_strandcorrectedsorted
open(READ, $inputfile) || die "Cannot open $inputfile!\n";
open(OUT, '>>', $outputfile) or die $!;
##read in line by line, check if the orthologs are longer than the phastcons and cut them in case
while(<READ>){
	chomp;

	if( /([^\s]+)\s+(\d+)\s+(\d+)\s+([^\s]+)\s+(\d+)\s+(\d+)\s+([+-])\s+([^\s]+)\s+([^\s]+)\s+(\d+)\s+(\d+)\s+([ACTGNactgn]+)\s+([+-])\s+([^\s]+)\s+([\d+-]+)\s+(\d+)/ ) { 
		#print "ori:$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$12\t$13\t$14\t$15\t$16\n";
		$phastchr = $1;
		$phaststart = $2;
		$phastend = $3;
		$hgchr = $4;
		$hgstart = $5;
		$hgend = $6;
		$hgstrand = $7;
		$orthgenome = $8;
		$orthchr = $9;
		$orthstart = $10;
		$orthend = $11;
		$orthseq = $12;
		$orthstrand = $13;
		$score = $14;
		$length = $15;
		$reallength = $16;
			
		if ( ($hgstart < $phaststart) && ($hgend > $phastend) ){
			$difference = $phaststart - $hgstart;
			#print "$difference \n";
                        $hgstart = $hgstart + $difference;
                        $orthstart = $orthstart + $difference;
                        if ( $orthstrand eq "+" ){
                                $orthseq =~ s/^.{$difference}//s;
                        }elsif( $orthstrand = "-" ){
                                $orthseq =~ s/.{$difference}$//s;
                        }
                        $reallength = $reallength - $difference;
			$difference = $hgend - $phastend;
			#print "$difference \n";
                        $hgend = $hgend - $difference;
                        $orthend = $orthend - $difference;
                        if ( $orthstrand eq "+" ){
                                $orthseq =~ s/.{$difference}$//s;
                        }elsif( $orthstrand = "-" ){
                                $orthseq =~ s/^.{$difference}//s;
                        }
                        $reallength = $reallength - $difference;
		}elsif( $hgstart < $phaststart ){
			$difference = $phaststart - $hgstart;
			#print "$difference \n";
			$hgstart = $hgstart + $difference;
			$orthstart = $orthstart + $difference;
			if ( $orthstrand eq "+" ){
				$orthseq =~ s/^.{$difference}//s;
			}elsif( $orthstrand = "-" ){
				$orthseq =~ s/.{$difference}$//s;
			}
			$reallength = $reallength - $difference;
		}elsif( $hgend > $phastend ){
			$difference = $hgend - $phastend;
			#print "$difference \n";
			$hgend = $hgend - $difference;
			$orthend = $orthend - $difference;
			if ( $orthstrand eq "+" ){
                                $orthseq =~ s/.{$difference}$//s;
                        }elsif( $orthstrand = "-" ){
                                $orthseq =~ s/^.{$difference}//s;
                        }
                        $reallength = $reallength - $difference;
		}
		print OUT "$1\t$2\t$3\t$4\t$hgstart\t$hgend\t$7\t$8\t$9\t$orthstart\t$orthend\t$orthseq\t$13\t$14\t$15\t$reallength\n";
	}else{
		print "Wrong format of line!\n $_\n";
	}
}

close(READ);
close(OUT);


	
