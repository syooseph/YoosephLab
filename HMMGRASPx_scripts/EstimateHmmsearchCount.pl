#!/usr/bin/perl -w
use strict;

my $file = shift;	# expecting hmmsearch result file
my %recruit_count;

open my $IN, "<$file" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  if(/^Query\:/)  {
    my @decom = split /\s+/, $_;
    my $name = $decom[1];
    my $count = 0;
    while(<$IN>)  {
      chomp;
      ++ $count;
      last if(!/\S+/);
    }
    $count = $count - 4;
    $count = 0 if($count < 0);
    $recruit_count{$name} = $count;
  }
}
close $IN;

foreach(sort keys %recruit_count)  {
  print "$_	$recruit_count{$_}\n";
}
