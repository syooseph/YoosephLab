#!/usr/bin/perl -w
use strict;

my $in_seq_file = shift;
my $out_seq_file = shift;

my $read_len = shift;
my $coverage = shift;

my @contents;
open my $IN, "<$in_seq_file" or die "Cannot open file: $!\n";
while(<$IN>) {
  chomp;
  push @contents, $_;
}
close $IN;

open my $OUT, ">$out_seq_file" or die "Cannot create file: $!\n";
my $index_global;
for($index_global = 0; $index_global < scalar(@contents); $index_global += 2) {

  my $info = $contents[$index_global];
  my $seq = $contents[$index_global + 1];
  
  if(length($seq) < $read_len)  {
    next;
  }

  $info =~ s/^>//g;
  chomp $info;
  chomp $seq;

  my $i;
  my $num_reads = int ($coverage * length($seq) / $read_len);
  
  for($i = 0; $i < $num_reads; ++ $i) {
    my $p = int(rand(length($seq)));
    #if($p - $read_len >= 0)  {
    #  my $start = $p - $read_len;
    #  my $end = $start + $read_len - 1;
    #  my $r_seq = substr($seq, $p - $read_len, $read_len);
    #  print $OUT ">$info $start:$end\n$r_seq\n";
    #}
    if($p + $read_len - 1 < length($seq))  {
      my $start = $p;
      my $end = $p + $read_len - 1;
      my $r_seq = substr($seq, $p, $read_len);
      print $OUT ">$info\:\:$start-$end\:\:N$i\n$r_seq\n";
    }  else  {
      -- $i;  # no read is generated, so reset the counter
    }
  }
}
