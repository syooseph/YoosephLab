#!/usr/bin/perl -w
use strict;

my $hmm_len = "/usr/local/projdata/0599/projects/SPA/czhong/Works/OralMeta/Results/long_HMM.list";
my $folder = "/usr/local/projdata/0599/projects/SPA/czhong/Works/OralMeta/Results/TargetedAssembly/GRASPxp";

my %len_hash;
open my $LIN, "<$hmm_len" or die "Cannot open file: $!\n";
while(<$LIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  $len_hash{$decom[0]} = $decom[1];
}
close $LIN;

#foreach(keys %len_hash) {
#  print "$_	$len_hash{$_}\n";
#}
#exit;

for(<$folder/*H*>)  {
  my $dir = $_;
  #print "$dir\n";
  chdir "$_" or die "Cannot change directory 1: $!\n";
  #next;
  foreach(<*>)  {
    my $name = $_;
    next if(!-d $name);
    #print "$name\n";
    chdir "$name" or die "Cannot change directory 2: $!\n";
    my $p_len = 0;
    my $n_len = 0;
    # computes peptide length
    if(-e "$dir/$name/SFA-SPA/spa.fasta")  {
      system "perl ~/ScriptRepo/ComputeN50.pl --fasta=$dir/$name/SFA-SPA/spa.fasta --min=60 >tmpN50";
      open my $IN, "<tmpN50" or die "Cannot open file: $!\n";
      my $line = <$IN>;
      chomp $line;
      $line =~ /is (\d+)\.$/;
      $p_len = $1 / $len_hash{$name};
      close $IN;
    }
    # computes nucleotide length
    if(-e "$dir/$name/SPAdes/contigs.fasta")  {
      system "perl ~/ScriptRepo/ComputeN50.pl --fasta=$dir/$name/SPAdes/contigs.fasta --min=180 >tmpN50";
      open my $IN, "<tmpN50" or die "Cannot open file: $!\n";
      my $line = <$IN>;
      chomp $line;
      $line =~ /is (\d+)\.$/;
      $n_len = $1 / ($len_hash{$name} * 3);
      close $IN;
    }
    print "$dir	$name	$p_len	$n_len\n";
    chdir "../" or die "Cannot change directory 3: $!\n";
  }
  chdir "../" or die "Cannot change directory 4: $!\n";
}
