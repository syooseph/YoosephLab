#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';

# batch submit SALT jobs

my $pfam_list = shift;		# a list of Pfam IDs of interest, one per line
my $data = shift;		# the FASTA format sequence file

$pfam_list = abs_path($pfam_list);
$data = abs_path($data);

my $SALT_dir = "/usr/local/projdata/0599/projects/SPA/czhong/Works/GRASPxp/Results/Glyco_SALT";
my $pfam_dir = "/usr/local/projdata/0599/projects/SPA/czhong/Data/Pfam_HMM";
my $source_dir = "/home/czhong/Tools/SALT";

open my $IN, "<$pfam_list" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  my $pfam_id = $_;
  mkdir "$SALT_dir/$pfam_id" or die "Cannot create directory: $!\n";
  chdir "$SALT_dir/$pfam_id" or die "Cannot create directory: $!\n";
  system "cp $pfam_dir/$pfam_id.hmm ./";
  system "nohup /usr/bin/time -v $source_dir/SALT.sh -m $pfam_id.hmm -f $data -o ./$pfam_id.out -E 10 &";
  
}
close $IN;
