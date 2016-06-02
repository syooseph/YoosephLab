#!/usr/bin/perl -w
use strict;

# the script is used to get species names out of the taxa mapping file
# working on directory /usr/local/projdata/0599/projects/SPA/czhong/Works/OralMeta/Results/UMap_ASsembled_BLAST

my $file = shift;

#foreach(</usr/local/projdata/0599/projects/SPA/czhong/Works/OralMeta/Results/UMap_ASsembled_BLAST/taxa*.list>)  {
  open my $IN, "<$file" or die "Cannot open file: $!\n";
  while(<$IN>)  {
    chomp;
    my @decom = split /\t/, $_;
    $decom[0] =~ s/\s+/\_/g;
    print "$decom[0]\n" if (!($decom[0] =~ /unmapped/));
  }
  close $IN;
#}
