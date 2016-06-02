#!/usr/bin/perl -w
use strict;

my $gene_id_list = shift;
my $id_pfam_map = shift;
my $pfam_dir = "/usr/local/depot/projects/SPA/czhong/Data/Pfam_full";

my %id_hash;
my %pfam_hash;

# get the gene ids with interest
open my $IN, "<$gene_id_list" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  $id_hash{$_} = 1;
}
close $IN;

# get the pfams of interest
open my $MIN, "<$id_pfam_map" or die "Cannot open file: $!\n";
while(<$MIN>)  {
  chomp;
  my @decom = split /[:\s]/, $_;
  if(exists $id_hash{$decom[1]})  {
    $pfam_hash{$decom[3]} = 1;
  }
}
close $MIN;

#foreach(keys %pfam_hash)  {
#  print "$_\n";
#}
#die;
# print the information
chdir "$pfam_dir" or die "Cannot change directory: pfam folder does not exist!";
foreach(<*>)  {
  my $pf_name = $_;
  #print "$pf_name\n";
  $pf_name =~ /(PF\d+)/;
  if(exists $pfam_hash{$1})  {
    system "ls $pf_name | wc -l";
  }
}
