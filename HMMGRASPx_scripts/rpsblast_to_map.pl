#!/usr/bin/perl -w
use strict;

my $rps_list = shift;		# the list of smp files that were used to construct the database
my $rps_results = shift;	# the RPS-blast results
my $cutoff = shift;		# the E-value cutoff

# construct ID mapping file
my %id_map;
open my $LIN, "<$rps_list" or die "Cannot open file: $!\n";
while(<$LIN>)  {
  chomp;
  open my $IN, "<$_" or die "Cannot open file: $!\n";
  my $id;
  my $desc;
  while(<$IN>)  {
    if(/\s+tag\s+id\s+(\d+)/)  {
      $id = 'gnl|CDD|' . $1;
    }  elsif(/\s+title\s+\S+\s+(\S+)\,/)  {
      $id_map{$id} = $1;
      last;
    }
  }
  close $IN;
}
close $LIN;

# process the RPS-blast results
my %best_score;
my %best_map;
open my $RIN, "<$rps_results" or die "Cannot open file: $!\n";
while(<$RIN>)  {
  chomp;
  if(/^\#\s(\d+)\shits\sfound/)  {
    if($1 > 0)  {
      my $line = <$RIN>;
      my @decom = split /\s+/, $line;
      #if($decom[10] < $cutoff)  {
      #  print "XX	$decom[0]	CHECKED||$id_map{$decom[1]}\n";
      #}
      if($decom[10] <= $cutoff and (!exists $best_score{$decom[0]} or $best_score{$decom[0]} > $decom[10]))  {
        $best_score{$decom[0]} = $decom[10];
        $best_map{$decom[0]} = $id_map{$decom[1]};
      }
    }
  }
}
close $RIN;

foreach(keys %best_map)  {
  print "XX	$_	CHECKED||$best_map{$_}\n";
}
