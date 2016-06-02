#!/usr/bin/perl -w
use strict;

my $annot_folder = shift;	# the folder contains the .kff files
my $tax_ref_map = shift;	# the mapping between tax ID and refseq NC ID
my $read_file = shift;		# the set of reads in FASTA format
my $hmm_parsed = shift;		# the parsed results of hmmer3 against full-length sequences
my $out_file = shift;		# the output file

# reads in annotations
my %annot_hash;
foreach(<$annot_folder/*.kff>)  {
  open my $IN, "<$_" or die "Cannot open file: $!\n";
  while(<$IN>)  {
    chomp;
    # index the annotation by gene ID
    my $line = $_;
    $line =~ /^(\S+)/;
    $annot_hash{$1} = $line;
  }
  close $IN;
}
print "Gene Annotation Information Loaded...\n";


# read in the mapping file
open my $MIN, "<$tax_ref_map" or die "Cannot open file: $!\n";
my %id_map;
while(<$MIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  $id_map{$decom[1]} = $decom[0];
}
close $MIN;
print "Taxanomic ID <-> RefSeq Mapping Information Loaded..\n";


# process and hash the reads according to their locations
open my $RIN, "<$read_file" or die "Cannot open file: $!\n";
my %read_hash;
while(<$RIN>)  {
  chomp;
  my $line = $_;
  if(/^>\S+ref\|(NC\_\d+)\S+\|\_(\d+)\_(\d+)\S+\/\d+\_(\d+)\_(\d+)\_(\S+)/)  {
    if(!exists $id_map{$1})  {
      next;
    }
    my $tax_id = $id_map{$1};
    my $begin;
    my $end;
    if($6 eq "+")  {
      $begin = $2 + $4 - 1;
      $end = $2 + $5 - 1;
    }  else  {
      $begin = $3 - $5 + 1;
      $end = $3 - $4 + 1;
    }
    my @loc = ($begin, $end, $line);
    push @{$read_hash{$tax_id}}, \@loc;
    #print "$tax_id	$begin	$end	$line	$1   $2      $3      $4      $5      $6\n";
  }
}
close $RIN;
print "Individual Reads Information Loaded...\n";

# sort the reads for faster access
foreach(keys %read_hash)  {
    my @array = @{$read_hash{$_}};
    @{$read_hash{$_}} = sort {$a->[0] <=> $b->[0]} @array;
}
print "Individual Reads Sorted...\n";

#foreach(keys %read_hash)  {
#  print "!!! $_ !!!\n";
#  my @array = @{$read_hash{$_}};
#  my $i; my $j;
#  for($i = 0; $i < scalar @array; ++ $i)  {
#    for ($j = 0; $j < 3; ++ $j)  {
#      print "$array[$i][$j]     ";
#    }
#    print "\n";
#  }
#}

# get regions identified by hmmer3
my %pf_id_hash;
my %inc_hash;
open my $HIN, "<$hmm_parsed" or die "Cannot open file: $!\n";
while(<$HIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  $pf_id_hash{$decom[1]} = $decom[0];
  my @decom2 = split /\s+/, $annot_hash{$decom[2]};
  $decom[2] =~ /^(\S+)\:/;
  my $tax_id = $1;
  $decom2[4] =~ /(\d+)\.\.(\d+)/;
  my $num1 = $1; my $num2 = $2;		# defines the regions of the gene
  if($decom2[4] =~ /complement/) {
     my $p = $num2;
     $num2 = $p - ($decom[3] * 3);
     $num1 = $p - ($decom[4] * 3);
  }  else  {
     my $p = $num1;
     $num1 = $p + ($decom[3] * 3);
     $num2 = $p + ($decom[4] * 3);
  }

  print "Region info:	$decom2[0]	$decom2[1]	$decom2[4]	$num1	$num2\n";
  
  # locate the position in the read list
  if(!exists($read_hash{$tax_id}))  {
    next;
  }
  my $left = 0;
  my $right = scalar(@{$read_hash{$tax_id}}) - 1;
  my $pivot;
  while($left < $right)  {
    $pivot = int(($left + $right) / 2);
    if(${$read_hash{$tax_id}}[$pivot][0] < $num1)  {
      $left = $pivot + 1;
    }  elsif(${$read_hash{$tax_id}}[$pivot][0] > $num1)  {
      $right = $pivot - 1;
    }  else  {
      last;
    }
  }
  my $search_begin;
  if($left < $right)  {
    $search_begin = $right - 1;
  }  else  {
    $search_begin = $pivot - 1;
  }
  if($search_begin < 0)  {
    $search_begin = 0;
  }
  # begin with the search
  #print "$tax_id	$num1	$search_begin	$read_hash{$tax_id}->[$search_begin][0]\n";
  for(my $k = $search_begin; $k < scalar(@{$read_hash{$tax_id}}); ++ $k)  {
    if(!($read_hash{$tax_id}->[$k][0] > $num2 || $read_hash{$tax_id}->[$k][1] < $num1))  {
      $inc_hash{$decom[1]}{$read_hash{$tax_id}->[$k][2]} = 1;
    }
    if($read_hash{$tax_id}->[$k][1] > $num2)  {
      last;
    }
  }
}
close $HIN;

open my $OUT, ">$out_file" or die "Cannot create file: $!\n";
foreach(keys %inc_hash)  {
  my $pfkey = $_;
  foreach(keys %{$inc_hash{$pfkey}})  {
    print $OUT "$pfkey	$pf_id_hash{$pfkey}	$_\n";
  }
}
close $OUT;

#foreach(keys %annot_hash)  {
#  print "$_	$annot_hash{$_}\n";
#}
