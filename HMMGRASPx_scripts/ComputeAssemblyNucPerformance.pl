#!/usr/bin/perl -w
use strict;

my $folder = "/usr/local/projdata/0599/projects/SPA/czhong/Works/OralMeta/Results/TargetedAssembly/HMMER3";
my $min_len = 180;
my $out = shift;

open my $OUT, ">$out" or die "Cannot create file: $!\n";

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
    if(-e "$dir/$name/SPAdes/contigs.fasta" && -s "$dir/$name/SPAdes/contigs.fasta")  {
      #system "perl ~/Tools/FragGeneScan1.20/run_FragGeneScan.pl --genome=SPAdes/contigs.fasta --out=spades_contig --complete=1 --train=complete";
      #system "~/Tools/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmsearch /usr/local/projdata/0599/projects/SPA/czhong/Works/OralMeta/Data/antiSMASH_ALL_HMM/$name.hmm spades_contig.faa >nuc_hmmer";
      system "perl ~/ScriptRepo/SummarizeHmmerResults.pl --results=nuc_hmmer --out=nuc_parsed";
      system "perl ~/ScriptRepo/FastaNewlineRemover.pl --fasta=$dir/$name/SPAdes/contigs.fasta --out=spades_contigs_no_newline.fasta";
      
      my %long_ids;
      open my $IN, "<$dir/$name/spades_contigs_no_newline.fasta" or die "Cannot open file: $!\n";
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
      open $IN, "<nuc_parsed" or die "Cannot open file: $!\n";
      while(<$IN>)  {
        chomp;
        my @decom = split /\s+/, $_;
        $decom[1] =~ /^(.*ID_\d+)\_/;
        ++ $num_true if ($decom[3] <= 0.01 && exists $long_ids{$1});
      }
      close $IN;
    }
    print $OUT "$dir	$name	$num_true	$total\n";
    chdir "../" or die "Cannot change directory 3: $!\n";
  }
  chdir "../" or die "Cannot change directory 4: $!\n";
}

close $OUT;
