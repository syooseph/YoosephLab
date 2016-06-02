#!/usr/bin/perl -w
use strict;

my $gene_sheet = shift;	# the list containing all genes of interest
my $raw_count_folder = shift;  # the folder containing raw counts of TP and total predictions

sub ComputeSenSpeFromArray($)  {
  # the input is expected to be a two dimensional array 
  # with at least one row, each row has at least 4 columns
  # each column has the following format:
  # ID true-size TP-tool1 Predicted-tool1 TP-tool2 Predicted-tool2 ...
  my $array = shift;
  my @performance;
  # initialize array, note that the first two columns do not correpond to performance
  for(my $i = 0; $i < scalar(@{$array->[0]}) - 2; ++ $i)  {
    $performance[$i] = 0;
  }
  # compute performance for each row
  my $total_entries = 0;
  for(my $i = 0; $i < scalar(@{$array}); ++ $i)  {
    for(my $j = 2; $j < scalar(@{$array->[$i]}); $j += 2)  {
      my $sen;
      if($array->[$i][1] == 0)  {
        $sen = 1;	# no true hit exist, so sensitivity is 1
      }  else  {
        $sen = $array->[$i][$j] / $array->[$i][1];
      }
      my $spe;
      if($array->[$i][$j + 1] == 0)  {
        $spe = 1;       # no hit identified, so specificity is 1
      }  else  {
        $spe = $array->[$i][$j] / $array->[$i][$j + 1];
      }
      # sum over existing performances
      $performance[$j - 2] += $sen;
      $performance[$j - 1] += $spe;
    }
    ++ $total_entries;
  }
  if($total_entries <= 0) {die "Empty matrix being input for performance computation!!!\n";}
  # compute average performance
  for(my $i = 0; $i < scalar(@performance); ++ $i)  {
    $performance[$i] /= $total_entries;
  }
  return @performance;
}

sub ComputeAUCFromMatrix($)  {
  # the input is a expected to be a two-dimensional matrix
  # with a format of:
  # sen_tool1 spe_tool1 sen_tool2 spe_tool2 ...
  # each row is expected to reflect the performance from E-value cutoff 10^-10 to 10
  # AUC computer based on simple extrapolation
  my $array = shift;
  my @aucs;
  my $num_rows = scalar(@{$array});
  for(my $j = 0; $j < scalar(@{$array->[0]}); $j += 2)  {
    # the first chunk is the triangle formed between the first data point and the origin
    my $area = 0.5 * $array->[0][$j] * (1 - $array->[0][$j + 1]);
    for(my $i = 1; $i < $num_rows; ++ $i)  {
      if($array->[$i][$j + 1] < $array->[$i - 1][$j + 1])  {
        # area is computable only when FDR is increasing
        $area += 0.5 * ($array->[$i - 1][$j + 1] - $array->[$i][$j + 1]) 
            * ($array->[$i - 1][$j] + $array->[$i][$j]);
      }
    }
    # compute the last area by extrapolation
    my $slope = ($array->[$num_rows - 1][$j] - $array->[$num_rows - 2][$j]) 
      / ($array->[$num_rows - 2][$j + 1] - $array->[$num_rows - 1][$j + 1]); 
    $slope /= 2; # extrapolated by taking the final slope 0
    my $intersect = $array->[$num_rows - 1][$j + 1] * $slope + $array->[$num_rows - 1][$j];
    $area += 0.5 * $array->[$num_rows - 1][$j + 1] * ($intersect + $array->[$num_rows - 1][$j]);
    push @aucs, $area;
  }
  return @aucs;
}

sub ComputeOptFScoreFromMatrix($)  {
  # the input is expected the same as ComputeAUCFromMatrix
  my $array = shift;
  my @f_scores;
  my $num_rows = scalar(@{$array});
  for(my $j = 0; $j < scalar(@{$array->[0]}); $j += 2)  {
    my $max_f = 0;
    for(my $i = 0; $i < $num_rows; ++ $i)  {
      my $f_score = 2 * $array->[$i][$j] * $array->[$i][$j + 1] / ($array->[$i][$j] + $array->[$i][$j + 1]);
      if($f_score > $max_f)  {
        $max_f = $f_score;
      }
    }
    push @f_scores, $max_f;
  } 
  return @f_scores;
}

# reads in the data sheet and define genes of interest
my %gid;
open my $IN, "<$gene_sheet" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  $gid{$_} = 1;
}
close $IN;

# read all raw count files from the "raw_count_folder"
# assuming E-value cutoffs range from 1e-10 to 1e1
my $e_min = -10; my $e_max = 1;
my @performances;
for(my $i = $e_min; $i <= $e_max; ++ $i)  {
  # assuming the files are named, i.e. performance.exp-5.count, in the folder
  open my $IN, "<$raw_count_folder/performance.exp$i.count" or die "Cannot open file: $!\n";
  my @count_matrix;
  while(<$IN>)  {
    chomp;
    my @decom = split /\s+/, $_;
    if(exists $gid{$decom[0]})  {	# only include the gene if it is in the gene sheet
      push @count_matrix, [@decom];
    }
  }
  close $IN;
  # compute the performance of the file
  my @perf = ComputeSenSpeFromArray(\@count_matrix);
  push @performances, [@perf];
}

my @fmeasure = ComputeOptFScoreFromMatrix(\@performances);
print "@fmeasure\n";
