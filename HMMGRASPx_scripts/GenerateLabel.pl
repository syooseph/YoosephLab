#!/usr/bin/perl -w
use strict;

my $file = shift;
my $res_label = "/usr/local/projdata/0599/projects/SPA/czhong/Works/OralRES/Data/resfams_merged.txt";

open my $IN, "<$res_label" or die "Cannot open file: $!\n";
my %hash;
while(<$IN>)  {
  chomp;
  my @decom = split /\t/, $_;
  $decom[5] =~ s/\s+/\_/g;
  $hash{$decom[2]} = $decom[5];
}
close $IN;

open my $FIN, "<$file" or die "Cannot open file: $!\n";
while(<$FIN>)  {
  chomp;
  my @decom = split /\s+/, $_;
  if(/^\#/)  {
    print "$decom[0]	$decom[1]	$decom[2]	$decom[3]	$decom[4]	$decom[5]	$decom[6]	$decom[7]	$decom[8]	$decom[9]	$decom[10]	$decom[11]	$decom[12]\n"
  }  else  {
    $decom[0] = $hash{$decom[0]} . '||' . $decom[0];
    print "$decom[0]	$decom[1]	$decom[2]	$decom[3]	$decom[4]	$decom[5]	$decom[6]	$decom[7]	$decom[8]	$decom[9]	$decom[10]	$decom[11]	$decom[12]\n";
  }
}
close $FIN;
