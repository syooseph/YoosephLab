#!/usr/bin/perl -w
use strict;

my $hmm_len = "/usr/local/projdata/0599/projects/SPA/czhong/Works/OralMeta/Results/long_HMM.list";
my $folder = "/usr/local/projdata/0599/projects/SPA/czhong/Works/OralMeta/Results/TargetedAssembly/HMMER3";
my $min_len = 60;

my %len_hash;
open my $LIN, "<$hmm_len" or die "Cannot open file: $!\n";
while(<$LIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  $len_hash{$decom[0]} = $decom[1];
}
close $LIN;


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
    my $num_true = 0;
    my $total = 0;
    # computes peptide length
    if(-e "$dir/$name/SFA-SPA/spa.fasta" && -s "$dir/$name/SFA-SPA/spa.fasta")  {
      #system "~/Tools/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmsearch /usr/local/projdata/0599/projects/SPA/czhong/Works/OralMeta/Data/antiSMASH_ALL_HMM/$name.hmm $dir/$name/SFA-SPA/spa.fasta >pep_hmmer";
      system "perl ~/ScriptRepo/SummarizeHmmerResults.pl --results=pep_hmmer --out=pep_parsed --cutoff=0.01";
      
      my %long_ids;
      open my $IN, "<$dir/$name/SFA-SPA/spa.fasta" or die "Cannot open file: $!\n";
      while(<$IN>)  {
        chomp;
        if(/^>(\S+)/)  {
          my $id = $1;
          my $next_seq = <$IN>;
          if (length($next_seq) - 1 >= $min_len)  {
            $long_ids{$id} = length($next_seq) - 1;
            ++ $total;
          }
        }
      }
      close $IN;
      open $IN, "<pep_parsed" or die "Cannot open file: $!\n";
      while(<$IN>)  {
        chomp;
        my @decom = split /\s+/, $_;
        ++ $num_true if ($decom[3] <= 10 && exists $long_ids{$decom[1]});
      }
      close $IN;
    }
    print "$dir	$name	$num_true	$total\n";
    chdir "../" or die "Cannot change directory 3: $!\n";
  }
  chdir "../" or die "Cannot change directory 4: $!\n";
}
