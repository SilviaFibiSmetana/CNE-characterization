#! /usr/bin/perl

use strict;
use warnings;

##usage: perl my_maftobed.pl maffile contained_genomes outputfile

##get files from commandline
unless($#ARGV == 2){
	die "usage: perl $0 maffile list_of_contained_genomes outputfile!\n";
}

my $maffile = $ARGV[0];
my $genomeslist = $ARGV[1];
my $outputfile = $ARGV[2];

##define variables 
#for the genomelist
my $numofgenomes = 470;
my @line;
my %genomeposition;
my $genome;
my $position;
#for the maf file
my $positioninvector;
my $chr = 0;
my $start;
my $end;
my $strand;
my $geningen;
my $chringen;
my $startingen;
my $seqingen;
my $strandingen;
my @geningen = ('?') x 470; 
my @chringen = ('?') x 470;
my @startingen = ('?') x 470;
my @seqingen = ('?') x 470;
my @strandingen = ('?') x 470;
my $score;
my @newline;


##read genomelist into hash
open(READ, $genomeslist) || die "Cannot open $genomeslist!\n";
while(<READ>){
	chomp;

	@line = split;
	#print "@line\n";
	$genome = $line[1];
	$position = $line[0];
	$genomeposition{$genome} = $position
}
close(READ);

#to check the hash
while(my($key, $value) = each %genomeposition){
	#print "$key, $value\n";
}


#read and process maffile
if ($maffile =~ /.gz$/) {
	open(READ, "gunzip -c $maffile|") || die "Cannot open pipe to $maffile!\n";
}
else {
	open(READ, $maffile) || die "Cannot open $maffile!\n";
}

while(<READ>){
	chomp;
	#print "line $_ \n";
	if( /^#/ ){
		next;
	}
	if( /^a score=/ ){
		@line = split /[\s=]+/;
		$score = $line[2];
		#print "here @line\n";
	}
	elsif (/^s\s+([^.]+)\.([^\s]+)\s+(\d+)\s+(\d+)\s+([+-])\s+(\d+)\s+([ACTGNactgn-]+)/){
		if( $1 eq 'hg38' ){
			#print "hgline @line \n";
			$chr = $2;
			$start = $3;
			$end = $3+$4-1; #is the -1 correct here? I would say yes, maftobed does not substract the 1.
			$strand = $5;
		}
		#print "there \n";
		$genome = $1;
		$positioninvector = $genomeposition{$genome};
	
		$geningen[$positioninvector] = $genome;
		$chringen[$positioninvector] = $2;
		$startingen[$positioninvector] = $3;
		$strandingen[$positioninvector] = $5;
		$seqingen[$positioninvector] = $7;
		#print "$positioninvector \n";
	}
	elsif( /^\s*$/ ){
		if ($chr eq '0'){
			#print "line1@line \n";
			#print "chr1 $chr \n";
			die "Wrong format!\n";
		}
		$geningen = join(',', @geningen);
		$chringen = join(',', @chringen);
		$startingen = join(',', @startingen);
		$seqingen = join(',', @seqingen);
		$strandingen = join(',', @strandingen);
		@newline = ($chr, $start, $end, $strand, $geningen, $chringen, $startingen, $seqingen, $strandingen, $score);
		open(OUT, '>>', $outputfile) or die $!;
		print OUT "@newline\n";
		close(OUT);
		#print "@newline\n";
		$chr = 0;
		$start = 0;
		$end = 0;
		$strand = 0;
		$genome = 0;
		$geningen = 0;
		$chringen = 0;
		$startingen = 0;
		$seqingen = 0;
		$strandingen = 0;
		$positioninvector = \undef;
		@geningen = ('?') x 470;
		@chringen = ('?') x 470;
		@startingen = ('?') x 470;
		@seqingen = ('?') x 470;
		@strandingen = ('?') x 470;
		$score = 0;
		@newline = 0;
	}
}
close(READ);

exit 0;
