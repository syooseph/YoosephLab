#!/usr/bin/perl -w
use strict;

my $in_directory = shift;
my $out_seq_file = shift;
my $read_len = shift;

if(-e "$out_seq_file")  {
  system "rm $out_seq_file";
}

foreach(<$in_directory/*.fa>) {
  system "perl RandomSlice.pl $_ temp $read_len";
  system "cat temp >>$out_seq_file"
}
