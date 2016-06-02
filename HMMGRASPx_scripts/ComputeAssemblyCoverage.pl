#!/usr/bin/perl -w
use strict;

my $hmm_len = "/usr/local/projdata/0599/projects/SPA/czhong/Works/OralMeta/Results/long_HMM.list";
my $folder = "/usr/local/projdata/0599/projects/SPA/czhong/Works/OralMeta/Results/TargetedAssembly/HMMER3";

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
    my $num_reads_p = 0;
    my $num_reads_n = 0;
    # computes peptide length
    if(-e "$dir/$name/SFA-SPA/spa.fasta")  {
      open my $IN, "<$dir/$name/SFA-SPA/spa.fasta" or die "Cannot open file: $!\n";
      while(<$IN>)  {
        chomp;
        if(/^>.*reads\:(\d+)/)  {
          my $nr = $1;
          my $tmp = <$IN>;
          $num_reads_p += $nr if (length($tmp) - 1 >= 60);
        }
      }
      close $IN;
    }
    # computes nucleotide length
    if(-e "$dir/$name/SPAdes/contigs.fasta")  {
      system "perl ~/ScriptRepo/FastaNewlineRemover.pl --fasta=$dir/$name/SPAdes/contigs.fasta --out=tmp_num_reads";
      open my $IN, "<tmp_num_reads" or die "Cannot open file: $!\n";
      while(<$IN>)  {
        chomp;
        if(/^>.*cov\_(.*)\_ID/)  {
          my $nr = $1;
          my $tmp = <$IN>;
          $num_reads_n += ($nr * (length($tmp) - 1) / 33) if (length($tmp) - 1 >= 180);
        }
      }
      close $IN;
    }
    print "$dir	$name	$num_reads_p	$num_reads_n\n";
    chdir "../" or die "Cannot change directory 3: $!\n";
  }
  chdir "../" or die "Cannot change directory 4: $!\n";
}
