#! /usr/bin/perl

use strict;
use warnings;

unless($#ARGV == 0){
	die "usage: perl $0 chr1short_intersected\n";
}

my $file = $ARGV[0]; 

my @line;
my @species;
my @chr;
my @start;
my @seq;
my @strand;
my $i;

open(READ, $file) || die "Cannot open $file!\n";
while(<READ>){
	chomp;

	@line = split;
	unless ( $#line == 13 ){
		die "$_: wrong format!\n";
	}

	@species = split(",", $line[7]);
	@chr = split(",", $line[8]);
	@start = split(",", $line[9]);
	@seq = split(",", $line[10]);
	@strand = split(",", $line[11]);

	for($i=0; $i<=$#species; $i++){
		if ($chr[$i] !~ m/^\?/){
			print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$species[$i]\t$chr[$i]\t$start[$i]\t$seq[$i]\t$strand[$i]\t$line[12]\t$line[13]\n";
		}
	}
}
close(READ);
