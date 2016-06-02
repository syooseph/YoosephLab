#!/usr/bin/perl -w
use strict;

my $blast_file = shift;	 # the full-sequence search results for the specific query sequence
#my $full_folder = shift; # the folder contains the blast full-sequence results
my $read_file = shift;   # the set of simulated read files
my $read_map_file = shift; # the mapping of the file to the orignal genome 
my $ev_cutoff = shift;	# evalue cutoff
my $out_dir = shift;	# the output directory

# reading information from the result that takes the sequence as query
my %aligned_location;
open my $IN, "<$blast_file" or die "Cannot open blast result file: $!\n";
while(<$IN>) {
  chomp;
  if(/^\#/)  {
    next;
  }
  my @decom = split /\s+/, $_;
  if(scalar(@decom) < 12 || $decom[11] > $ev_cutoff)  {
    next;
  }
  push @{$aligned_location{$decom[0]}}, [$decom[1], $decom[8], $decom[9] - 1, $decom[10] - 1, $decom[11]];
  #if(exists $aligned_location{$decom[0]})  {
    #push @{$aligned_location{$decom[0]}}, $decom[1];	# the target genome name
    #push @{$aligned_location{$decom[0]}}, $decom[8];	# the frame
    #push @{$aligned_location{$decom[0]}}, $decom[9] - 1;	# the begin
    #push @{$aligned_location{$decom[0]}}, $decom[10] - 1;	#the end
    #push @{$aligned_location{$decom[0]}}, $decom[11];	# the evalue
    
  #} else  {
  #  $aligned_location{$decom[0]}= [{$decom[1]}, $decom[8], $decom[9] - 1, $decom[10] - 1, $decom[11]];
  #}
}
close $IN;

#foreach(keys %aligned_location)  {
#  print "*******************$_\n";
#  for(my $k = 0; $k < scalar(@{$aligned_location{$_}}); ++ $k) {
#    print "@{$aligned_location{$_}[$k]}\n";
#  }
#}
#die;
# reading information from the original read file
my $index = 0;
my %read_id_hash;
open my $RIN, "<$read_file" or die "Cannot open file: $!\n";
while(<$RIN>)  {
  chomp;
  if(/^>(.*)/)  {
    $read_id_hash{$1} = $index ++;
  }
}
close $RIN;

# reading where each read is mapped to
open my $MIN, "<$read_map_file" or die "Cannot open file: $!\n";
my %read_map_info;
while(<$MIN>) {
  chomp;
  my @decom = split /\s+/, $_;
  push @{$read_map_info{$decom[1]}}, [@decom];
}
close $MIN;
#die;


#foreach(keys %aligned_location) {
#  print "$_\n";
#  print "@{$aligned_location{$_}}\n";
#}
#die;
my $i;
my $j;
foreach(keys %aligned_location)  {
  my $num_reads = 0;
  my $query_id = $_;
  #my @included_reads;
  print "Working on $query_id......\n";
  open my $OUT, ">$out_dir/$query_id.ref" or die "Cannot create file: $!\n";
  for($i = 0; $i < scalar(@{$aligned_location{$query_id}}); ++ $i)  {
    my $genome_id = $aligned_location{$query_id}[$i][0];
    for($j = 0; $j < scalar(@{$read_map_info{$genome_id}}); ++ $j) {
      if( (($read_map_info{$genome_id}[$j][2] eq '-' && $aligned_location{$query_id}[$i][1] < 0) || 
          ($read_map_info{$genome_id}[$j][2] eq '+' && $aligned_location{$query_id}[$i][1] > 0))) {
        #print "$read_map_info[$j][3]	$read_map_info[$j][4]	$aligned_location{$query_id}[$i][2]	$aligned_location{$query_id}[$i][3]\n";
        if($read_map_info{$genome_id}[$j][2] eq '+' &&
           $read_map_info{$genome_id}[$j][3] >= $aligned_location{$query_id}[$i][2] && 
           $read_map_info{$genome_id}[$j][4] <= $aligned_location{$query_id}[$i][3])  {
          
          # the read is included
          #push @{$included_reads{$query_id}}, $read_map_info[$j][0];
          print $OUT "$read_map_info{$genome_id}[$j][0]	$read_id_hash{$read_map_info{$genome_id}[$j][0]}\n";
        } elsif($read_map_info{$genome_id}[$j][2] eq '-' &&
           $read_map_info{$genome_id}[$j][3] <= $aligned_location{$query_id}[$i][2] &&
           $read_map_info{$genome_id}[$j][4] >= $aligned_location{$query_id}[$i][3])  {
	  print $OUT "$read_map_info{$genome_id}[$j][0]     $read_id_hash{$read_map_info{$genome_id}[$j][0]}\n";
        }
      }
    }
  }
  close $OUT;
}
