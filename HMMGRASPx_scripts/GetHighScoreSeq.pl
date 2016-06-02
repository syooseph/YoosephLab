#!/usr/bin/perl -w
use strict;

my $folder = "/usr/local/projdata/0599/projects/SPA/czhong/Works/OralMeta/Results/TargetedAssembly/GRASPxp";
my $out = shift;

#open my $OUT, ">$out" or die "Cannot create file: $!\n";

for(<$folder/*H*>)  {
  my $dir = $_;
  #print "$dir\n";
  chdir "$_" or die "Cannot change directory 1: $!\n";
  $dir =~ /(\dH\d)$/;
  my $time = $1;
  #next;
  foreach(<*>)  {
    my $name = $_;
    next if(!-d $name);
    #print "$name\n";
    chdir "$name" or die "Cannot change directory 2: $!\n";
    my $num_true = 0;
    my $total = 0;
    # computes peptide length
    #print "$dir/$name/SPAdes/contigs.fasta\n";
    if(-e "$dir/$name/SPAdes/contigs.fasta" && -s "$dir/$name/SPAdes/contigs.fasta")  {
      #system "perl ~/Tools/FragGeneScan1.20/run_FragGeneScan.pl --genome=SPAdes/contigs.fasta --out=spades_contig --complete=1 --train=complete";
      #system "~/Tools/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmsearch /usr/local/projdata/0599/projects/SPA/czhong/Works/OralMeta/Data/antiSMASH_ALL_HMM/$name.hmm spades_contig.faa >nuc_hmmer";
      system "perl ~/ScriptRepo/SummarizeHmmerResults.pl --results=nuc_hmmer --out=nuc_parsed";
      #print "$dir	$name\n";
      system "cat nuc_parsed \| awk \'\{if\(\$4<=0.01\) print \$2\}\' >id_list";
      #system "cat id_list";
      system "perl ~/ScriptRepo/GetFastaSeqByID.pl --fasta=spades_contig.ffn --id=id_list --append=$time\_$name --out=fetched.seq";
      system "cat fetched.seq >>$out";
    }
#    print $OUT "$dir	$name	$num_true	$total\n";
    chdir "../" or die "Cannot change directory 3: $!\n";
  }
  chdir "../" or die "Cannot change directory 4: $!\n";
}

#close $OUT;
