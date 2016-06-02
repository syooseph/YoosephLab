#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;

foreach my $f (@ARGV) {
    my $in = Bio::SeqIO->newFh(-file=>$f, -format=>"genbank");
    while (my $seq = <$in>) {
        my $source_acc=$seq->id;
        foreach my $feat ($seq->get_SeqFeatures() ) {
            next if ($feat->primary_tag ne 'CDS');
            my $acc;
            if ($feat->has_tag('locus_tag')) {
                ($acc) = $feat->get_tag_values('locus_tag');
            }
            else {
                ($acc) = $feat->get_tag_values('protein_id');
            }
            my ($description) = $feat->get_tag_values('product');
            my ($seq) = $feat->get_tag_values('translation');
            print ">$acc $description\n$seq\n";
        }
    }
}
