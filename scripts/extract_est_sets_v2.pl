#! /usr/bin/perl

use strict;
use warnings;

##usage: perl extract_est_sets.pl CLADEESTNUMBERS ALLORTHOLOGSNUMBERED CLADEESTSETS

##get files from commandline
unless($#ARGV == 2){
	die "usage: perl $0 CLADEESTNUMBERS ALLORTHOLOGSNUMBERED CLADEESTSETS \n";
}

my $cladenumfile = $ARGV[0];
my $inputfile = $ARGV[1];
my $outputfile = $ARGV[2];

#warn "# Assuming that $inputfile is sorted, so that elements of a set occur in consecutive lines. Also, assuming one line per assembly.\n";

my $assembly1 = 'hg38';
my @array;
my $line;

my %numhash;

##-------------------------------------------------
# READ CLADE ASSIGNMENTS INTO an ARRAY
# AND TRANSFER INTO HASH FOR EASY CHECKUP
##-------------------------------------------------

open (READ, $cladenumfile) || die "Cannot open $cladenumfile!\n"; #my $fh, '<', 'input.txt' or die "Cannot open file: $!";
while (my $line = <READ>) {
    chomp($line);  # Remove trailing newline
    push @array, $line;  # Add the line to the array
}
close READ;
#print "File contents: @array\n";

%numhash = map { $_ => 1 } @array;

##-------------------------------------------------
# SELECT THE HUMAN GENOME LINES THAT ARE THE SETS
# THAT ARE IN THE ESTIMATED LIST
##-------------------------------------------------

# OPEN OUTPUT FILE
open(OUT, '>', $outputfile) or die $!;

# OPEN INPUT FILE
open(READ, $inputfile) || die "Cannot open $inputfile!\n";
while(<READ>){
	chomp;
	##process all the other lines
	if( /^([^\s]+)\s+(\d+)\s+(\d+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+(\d+)\s+(\d+)\s+([^\s]+)\s+([^\s]+)\s+(\d+)\s+([^\s]+)\s+(\d+)/ ) {
		#print "$7\t$15\n";
		#If the line is a human line and is part of the estimated sets, print it
		if(($7 eq $assembly1) && (exists $numhash{$15})) {
			#print "$7\t$15\n";
			print OUT "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$12\t$13\t$14\t$15\n";
		}else{
			##do not do anything
			#print "not in set\n";
		}
	}else{	
		die "$_: Wrong format of line!\n";
	}
}
close(READ);
close(OUT);

exit 0;
