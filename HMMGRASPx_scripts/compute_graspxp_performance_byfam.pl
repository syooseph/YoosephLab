#!/usr/bin/perl -w
use strict;

my $tp_file = shift;		# the annotation of true-positive reads
my $read_file = shift;		# the fasta file containing the original reads
my $graspxp_mapping = shift;	# the mapping file generated by GRASPxp

# process the tp_file
my %ref_hash;
open my $TIN, "<$tp_file" or die "Cannot open file: $!\n";
while(<$TIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  $ref_hash{$decom[1]}{$decom[2]} = 0;
}
close $TIN;

# process the fasta file 
my %read_hash;
open my $RIN, "<$read_file" or die "Cannot open file: $!\n";
while(<$RIN>)  {
  chomp;
  if(/^>/)  {
    $read_hash{$_} = 0;
  }
}
close $RIN;

# process the grasp mapping data
my %tp;
my %fp;
my %fn;
open my $IN, "<$graspxp_mapping" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  my $key = '>' . $decom[1];
  $decom[2] =~ /\|\|(\S*)/;
  my $fam_desc = $1;
  # check if the read is correctly assigned
  if(exists $ref_hash{$fam_desc}{$key})  {
    # TP
    $tp{$fam_desc}{$key} = 1;
    $ref_hash{$fam_desc}{$key} = 1;
  }  else  {
    # FP
    $fp{$fam_desc}{$key} = 1;
  }
  $read_hash{$key} = 1;
}
close $IN;

foreach(keys %ref_hash)  {
  my $fam = $_;
  foreach(keys %{$ref_hash{$fam}})  {
    if($ref_hash{$fam}{$_} == 0)  {
      # FN
      $fn{$fam}{$_} = 1;
    }
    $read_hash{$_} = 1;
  }
}

foreach(keys %read_hash)  {
  if($read_hash{$_} == 1)  {
    # remove non-TNs
    delete $read_hash{$_};
  }
}

foreach(keys %tp)  {
  my $numtp = scalar(keys %{$tp{$_}});
  my $numfp = scalar(keys %{$fp{$_}});
  my $numfn = scalar(keys %{$fn{$_}});

  print "$_	$numtp	$numfp	$numfn\n";
}
