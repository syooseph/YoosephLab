#!/usr/bin/perl -w
use strict;

my $genome_ids = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Data/marine.core.id";
my $mapping_file = "/usr/local/projects/DB/kegg/kegg_current/genes/genome/genome_refseq.list";
my $taxa_file = "/usr/local/projects/DB/kegg/kegg_current/genes/taxonomy";
my $base_folder = "/usr/local/projects/DB/kegg/kegg_current/genes/organisms";
my $output_folder = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Data/Marine_Core_Genomes/Annotation";

# load into the refseq genome ids
my %grefseq_id;
open my $IIN, "<$genome_ids" or die "Cannot open file: $!\n";
while(<$IIN>)  {
  chomp;
  $grefseq_id{$_} = 1;
}
close $IIN;

# find the taxa IDs
my %taxa_id;
open my $MIN, "<$mapping_file" or die "Cannot open file: $!\n";
while(<$MIN>)  {
  chomp;
  /^genome\:(.*)\s+rs\:(.*)/;
  if(exists $grefseq_id{$2})  {
    $taxa_id{$1} = $2;
  }
}
close $MIN;

#foreach(keys %taxa_id)  {
#  print "$_	$taxa_id{$_}\n";
#}

# refers to the taxa file and copy the protein sequences to output folder
open my $TIN, "<$taxa_file" or die "Cannot open file: $!\n";
while(<$TIN>)  {
  chomp;
  if(/^#/)  {
    next;
  }
  my @decom = split /\s+/, $_;
  if(exists $taxa_id{$decom[2]})  {
    print "Copying peptide and annotation for genome $decom[2]...\n";
    system "cp $base_folder/$decom[1]/$decom[2].pep $output_folder/$decom[2]-$taxa_id{$decom[2]}.pep";
    system "cp $base_folder/$decom[1]/$decom[2].kff $output_folder/$decom[2]-$taxa_id{$decom[2]}.kff";
  }
}
close $TIN;
