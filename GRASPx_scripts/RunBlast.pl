#!/usr/bin/perl -w
use strict;

my $query_seq_dir = shift;
my $db_name = shift;
my $out_dir = shift;

my $large_num = 1000000;

foreach(<$query_seq_dir/*.fa>) {
  $_ =~ /(.*)\/(.*)/;
  #system "blastp -query $_ -db $db_name -num_descriptions $large_num -num_alignments $large_num >$out_dir/$2.blast.alignment";
  system "blastp -query $_ -db $db_name -num_descriptions $large_num -num_alignments $large_num -outfmt 7 -num_threads 2 >$out_dir/$2.blast.summary";
}
