#!/usr/bin/perl -w
use strict;

my $file = shift;
my $out_folder = shift;

open my $IN, "<$file"  or die "Cannot open file: $!\n";

while(<$IN>)  {
  chomp;
  if(/^>/)  {
    my $header = $_;
    my @decom = split /\s+/, $_;
    $decom[0] =~ s/^>//g;
    $decom[0] =~ s/[\.\|]/\_/g;
    open my $OUT, ">$out_folder/$decom[0].fna" or die "Cannot create file: $!\n";
    print $OUT "$header\n";
    my $temp = <$IN>;
    chomp $temp;
    print $OUT "$temp";
    close $OUT;
  }
}
close $IN;
