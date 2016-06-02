#!/usr/bin/perl -w
use strict;

my $hmmer_file = shift;		# the hmmsearch output file
my $cutoff = shift;		# the E-value cutoff

my %best_score;
my %best_map;
open my $IN, "<$hmmer_file" or die "Cannot open file: $!\n";
my $PF_id;
my $desc;
my $read;
while(<$IN>)  {
  chomp;
  if(/^Accession\:\s+(\S+)/)  {
    $PF_id = $1;
    next;
  }  elsif(/^Query\:\s+(\S+)/)  {
    $desc = $1;
    next;
  }  elsif(/^\>\>\s+(\S+)/)  {
    $read = $1;
  }
  #my @decom = split /\s+/, $_;
  #if(scalar(@decom) == 10 and $decom[1] =~ /\d+/)  {
  #  if($decom[1] < $cutoff)  {
  #    print "XX	$decom[9]	CHECKED||$desc\n";
  #  }
  #}

  if(/\?/ || /\!/)  {
    my @decom = split /\s+/, $_;
    if($decom[6] <= $cutoff and (!exists $best_score{$read} || $best_score{$read} > $decom[6]))  {
      $best_score{$read} = $decom[6];
      $best_map{$read} = $desc;
    }
  }
}
close $IN;

foreach(keys %best_map)  {
  print "XX	$_	CHECKED||$best_map{$_}\n";
}
