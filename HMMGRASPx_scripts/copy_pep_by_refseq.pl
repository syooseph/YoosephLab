#!/usr/bin/perl -w
use strict;

my $file = shift;
open my $IN, "<$file" or die "Cannot open file: $!\n";

my %refseq_id;
while(<$IN>)  {
  chomp;
  my $key = 'rs:' . $_;
  $refseq_id{$key} = 1;
}
close $IN;

open my $GIN, "</usr/local/db/kegg/kegg_current/genes/genome/genome_refseq.list" or die "Cannot open file: $!\n";
my %taxa_id;
while(<$GIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  if(exists $refseq_id{$decom[1]})  {
    $decom[0] =~ /genome:(.*)/;
    $taxa_id{$1} = $decom[1];
  }
}
close $GIN;

open my $TIN, "</usr/local/db/kegg/kegg_current/genes/taxonomy" or die "Cannot open file: $!\n";
while(<$TIN>)  {
  chomp;
  if(/^\#/)  {
    next;
  }
  my @decom = split /\s+/, $_;
  if(exists $taxa_id{$decom[2]})  {
    print "$decom[1]  $taxa_id{$decom[2]}\n";
    system "cp /usr/local/db/kegg/kegg_current/genes/organisms/$decom[1]/*.pep /usr/local/db/kegg/kegg_current/genes/organisms/$decom[1]/*.kff /usr/local/db/kegg/kegg_current/genes/organisms/$decom[1]/*pfam.list /usr/local/projdata/0599/projects/SPA/czhong/Works/GRASPxp/Data/"
  }
}
