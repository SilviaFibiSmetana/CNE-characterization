use strict;
use warnings;


unless ($#ARGV == 3) {
    die "usage: perl $0 bed output\n";
}

my $targetdistance = $ARGV[0];
my $phastconfraction = $ARGV[1];
my $bed = $ARGV[2];
my $outputfile = $ARGV[3];

print STDOUT "# Input: $bed\n";

my @element_array;
my ($chr1, $st1, $end1);
my ($chr2, $st2, $end2);

# Function to process and merge elements
sub process_elements {
    my ($element_array_ref, $chr, $output_handle) = @_;
    my @element_array = @$element_array_ref;
    
    # CHECK IF THE ARRAY CONTAINS ONLY ONE PHASTCON AND IF JUMP OVER THE MERGING
    if (scalar @element_array == 1) {
    	#print "is one \n";
        my ($start, $end) = @{$element_array[0]};
        print $output_handle "$chr\t$start\t$end\n";
    }else{
        #print "is more than one \n";
    
    
    
        # SORT ELEMENTS BY COORDINATES
        @element_array = sort { $a->[0] <=> $b->[0] or $a->[1] <=> $b->[1] } @element_array;

        # CALCULATE DISTANCES BETWEEN ELEMENTS
        my (@distances, @elements_index);
        for (my $i = 1; $i <= $#element_array; $i++) {
            my $dist = $element_array[$i]->[0] - $element_array[$i - 1]->[1];
            die "Error. Elements should not overlap!\n" if $dist < 0;
            push @distances, $dist;
            push @elements_index, [$i - 1, $i];
        }

        # SORT DISTANCES
        my @distances_index = sort { $distances[$a] <=> $distances[$b] } 0 .. $#distances;
        @distances = @distances[@distances_index];
        @elements_index = @elements_index[@distances_index];

        # MERGE ELEMENTS ACCORDING TO DISTANCES
        my (%elements_ref, %elements_merged_to_ref, %phastCons_content);
        for (my $i = 0; $i <= $#distances; $i++) {
            my ($left, $right) = @{$elements_index[$i]};
            my $ref = defined $elements_ref{$left} ? $elements_ref{$left} : $left;

            # Calculate PhastCons content
            my $phastCons_content = ($phastCons_content{$ref} // $element_array[$ref]->[1] - $element_array[$ref]->[0] + 1) +
                                ($phastCons_content{$right} // $element_array[$right]->[1] - $element_array[$right]->[0] + 1);

            # Calculate merged length and ratio
            my $left_border  = $element_array[$ref]->[0];
            my $right_border = $element_array[$right]->[1];
            my $merged_length = $right_border - $left_border;
            my $phastCons_ratio = $phastCons_content / $merged_length;

            if ($phastCons_ratio > $phastconfraction) {
                $elements_ref{$right} = $ref;
                $elements_merged_to_ref{$ref}{$right}++;
                if (defined $elements_merged_to_ref{$right}) {
                    foreach my $merged (keys %{$elements_merged_to_ref{$right}}) {
                        $right_border = $element_array[$merged]->[1] if $element_array[$merged]->[1] > $right_border;
                        $elements_ref{$merged} = $ref;
                        $elements_merged_to_ref{$ref}{$merged}++;
                    }
                }
                $phastCons_content{$ref} = $phastCons_content;
            }
        }

        # Collect and print final merged elements
        my %final_elements;
        foreach my $i (0 .. $#element_array) {
            my $ref = defined $elements_ref{$i} ? $elements_ref{$i} : $i;
            my $left_border  = $element_array[$ref]->[0];
            my $right_border = $element_array[$ref]->[1];
            if (defined $elements_merged_to_ref{$ref}) {
                foreach my $merged (keys %{$elements_merged_to_ref{$ref}}) {
                    $right_border = $element_array[$merged]->[1] if $element_array[$merged]->[1] > $right_border;
                }
            }
            $final_elements{$ref} = [$left_border, $right_border];
        }

        foreach my $ref (sort { $a <=> $b } keys %final_elements) {
            my ($start, $end) = @{$final_elements{$ref}};
            print $output_handle "$chr\t$start\t$end\n";
        }
    }
}

#Main part

open(my $OUT, '>', $outputfile) or die $!;

my $first_line_done = 0;

open(READ, $bed) or die "Cannot open $bed!\n";
while (<READ>) {
    chomp;
    #Read the file line by line to compare two consecutive lines
    if (!$first_line_done) {
        if (/^([^\t]+)\t([^\t]+)\t([^\t]+)/) {
            ($chr1, $st1, $end1) = ($1, $2, $3);
            push @element_array, [$st1, $end1];
            $first_line_done = 1;
        } else {
            die "Wrong format of first line!\n";
        }
    } else {
        if (/^([^\t]+)\t([^\t]+)\t([^\t]+)/) {
            ($chr2, $st2, $end2) = ($1, $2, $3);
            #If the distance between two elements is smaller than 50 bp, they belong to one cluster to investigate with the function and they are collected in an array
            if (($st2 - $end1) <= $targetdistance && $chr2 eq $chr1) {
                push @element_array, [$st2, $end2];
                ($chr1, $st1, $end1) = ($chr2, $st2, $end2);
            } else { #if a cluster ends its processed with the function and the array is emptied and started new
                process_elements(\@element_array, $chr1, $OUT);
                @element_array = ([$st2, $end2]);
                ($chr1, $st1, $end1) = ($chr2, $st2, $end2);
                    
            }
            if (eof) { #at the very end of the file the last cluster is also processed with the function
                process_elements(\@element_array, $chr1, $OUT);
            }
        } else {
            die "Wrong format of line!\n";
        }
    }
}
close(READ);
close($OUT);

exit 0;
