#!/usr/bin/perl -w
use strict;

my $hmm_full_file = shift;		# the hmmsearch result file
my $peptide_file = shift;		# the target peptide file in multi-fasta format

open my $IN, "<$hmm_full_file" or die "Cannot open file: $!\n";

my $name;
my $id;
my $target;
my @identified_regions;
while(<$IN>)  {
  chomp;
  if(/^Query:\s+(\S+)/)  {
    $name = $1;
  } elsif(/^Accession:\s+(PF\d+)\./)  {
    $id = $1;
  } elsif(/^>>\s+(\S+)/)  {
    $target = $1;
  } elsif(/\!/)  {
    my @decom = split /\s+/, $_;
    my @info = ($name, $id, $target, $decom[10], $decom[11]);
    push @identified_regions, \@info;
    #print "$name, $id, $target, $decom[10], $decom[11]\n";
  }
}

my %seq_hash;
open my $PIN, "<$peptide_file" or die "Cannot open file: $!\n";
while(<$PIN>)  {
  chomp;
  if(/^>(\S+)/)  {
    my $id = $1;
    my $line = <$PIN>;
    chomp $line;
    $seq_hash{$id} = $line;
  }
}
close $PIN;

foreach(@identified_regions)  {
  my @info = @{$_};
  if(exists $seq_hash{$info[2]})  {
    my $seq = substr($seq_hash{$info[2]}, $info[3] - 1, $info[4] - $info[3] + 1);
    print ">contig||$info[2]||$info[0]||$info[1]\n$seq\n";
  }
}
