use strict;
use warnings;

unless ( $#ARGV == 2 ){
	die "usage: perl $0 list_of_assemblies path_to_assemblies alignments!\n";
}

my $list_of_assemblies = $ARGV[0];
my $path = $ARGV[1];
my $alignments = $ARGV[2];

my %chromInfo;
my ($assembly, $chr, $size);

%chromInfo = %{read_assemblies($list_of_assemblies, $path)};

my @line;
my $chrlength;
my $start;

open(READ, $alignments) || die "Cannot open $alignments!\n";
while(<READ>){
	chomp;

	@line = split;
	unless ( $#line == 13 ){
		die "$_: wrong format!\n";
	}

	if ($line[11] =~ "-"){
		$assembly = $line[7];
		$chr = $line[8];
		$chrlength = $chromInfo{$assembly}{$chr};
		$start = $chrlength - $line[9];
		print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$start\t$line[10]\t$line[11]\t$line[12]\t$line[13]\n";
	}else{
		print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[10]\t$line[11]\t$line[12]\t$line[13]\n";
	}
}
close(READ);

sub read_assemblies {
	my $file = $_[0];
	my $path = $_[1];

	my $assembly;
	my $assembly_file;
	my ($chr, $size);

	open(READ, $file) || die "Cannot open $file!\n";
	while(<READ>){
		chomp;
		$assembly = $_;
		$assembly_file = $path.$assembly."/chromInfo.txt";
		unless( -e $assembly_file ){
			next;
		}
		open(READ_ASSEMBLY, $assembly_file) || die "Cannot open $assembly_file!\n";
		while(<READ_ASSEMBLY>){
			chomp;
			if ( /^([^\s]+)\s(\d+)/ ){
				$chr = $1;
				$size = $2;
				$chromInfo{$assembly}{$chr} = $size;
			}
		}
		close(READ_ASSEMBLY);
	}
	close(READ);
	
	return(\%chromInfo);
}

