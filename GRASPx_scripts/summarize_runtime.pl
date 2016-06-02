#!/usr/bin/perl -w
use strict;

my $query_file = shift;
my $gx_log_file = shift;

my %len_hash;
open my $QIN, "<$query_file" or die "Cannot open file: $!\n";
while(<$QIN>)  {
  chomp;
  if(/^>/)  {
    />(.*)/;
    my $name = $1;
    my $temp = <$QIN>;
    $len_hash{$name} = length($temp);
  }
}
close $QIN;

my %runtime_hash;
open my $GIN, "<$gx_log_file" or die "Cannot open file: $!\n";
while(<$GIN>)  {
  chomp;
  if(/\[(.*)\].*\[(.*)\]/)  {
    #print "$1	$2\n";
    my $name = $1;
    my @decom = split /\:/, $2;
    # convert the time into secs
    my $time = $decom[0] * 3600 + $decom[1] * 60 + $decom[2];
    $runtime_hash{$name} += $time;
  }
}
close $GIN;

foreach(sort keys %runtime_hash)  {
  print "$_  $len_hash{$_}  $runtime_hash{$_}\n";
}
