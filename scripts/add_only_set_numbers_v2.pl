#! /usr/bin/perl

use strict;
use warnings;

##usage: perl add_only_set_numbers.pl inputfile  outputfile 

##get files from commandline
unless($#ARGV == 1){
	die "usage: perl $0 inputfile outputfile!\n";
}

my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];


my $assembly1 = 'hg38';
my $assembly2;
my $setnumber = 0;
my ($line, $chunk);

open(OUT, '>', $outputfile) or die $!;

open(READ, $inputfile) || die "Cannot open $inputfile!\n";


while(<READ>){
	chomp;
	##process the lines
	if( /^([^\s]+)\s+(\d+)\s+(\d+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+(\d+)\s+(\d+)\s+([^\s]+)\s+([^\s]+)\s+(\d+)\s+([^\s]+)/ ) {
		$assembly2 = $7;
		#print "$assembly2\n";
		$line = $_;
		#print "$line\n";				
		
		if($assembly2 eq "hg38") {
			$setnumber++;
		}else{
			#dont do anything
		}
		$chunk = $line."\t".$setnumber."\n";
		print OUT $chunk;	
	}else{	
		die "$_: Wrong format of line!\n";
	}
}

close(READ);
close(OUT);

exit 0;
