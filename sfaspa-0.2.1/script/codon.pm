#!/usr/bin/env perl

use strict;
use warnings;
use Set::Scalar;

package codon;

my %CodonTable = (
				  "TTT" => "F",
				  "TTC" => "F",
				  "TTA" => "L",
				  "TTG" => "L",
				  
				  "TCT" => "S",
				  "TCC" => "S",
				  "TCA" => "S",
				  "TCG" => "S",

				  "TAT" => "Y",
				  "TAC" => "Y",
				  "TAA" => "*",
				  "TAG" => "*",

				  "TGT" => "C",
				  "TGC" => "C",
				  "TGA" => "*",
				  "TGG" => "W",

				  "CTT" => "L",
				  "CTC" => "L",
				  "CTA" => "L",
				  "CTG" => "L",
				  
				  "CCT" => "P",
				  "CCC" => "P",
				  "CCA" => "P",
				  "CCG" => "P",

				  "CAT" => "H",
				  "CAC" => "H",
				  "CAA" => "Q",
				  "CAG" => "Q",

				  "CGT" => "R",
				  "CGC" => "R",
				  "CGA" => "R",
				  "CGG" => "R",

				  "ATT" => "I",
				  "ATC" => "I",
				  "ATA" => "I",
				  "ATG" => "M",
				  
				  "ACT" => "T",
				  "ACC" => "T",
				  "ACA" => "T",
				  "ACG" => "T",

				  "AAT" => "N",
				  "AAC" => "N",
				  "AAA" => "K",
				  "AAG" => "K",

				  "AGT" => "S",
				  "AGC" => "S",
				  "AGA" => "R",
				  "AGG" => "R",

				  "GTT" => "V",
				  "GTC" => "V",
				  "GTA" => "V",
				  "GTG" => "V",
				  
				  "GCT" => "A",
				  "GCC" => "A",
				  "GCA" => "A",
				  "GCG" => "A",

				  "GAT" => "D",
				  "GAC" => "D",
				  "GAA" => "E",
				  "GAG" => "E",

				  "GGT" => "G",
				  "GGC" => "G",
				  "GGA" => "G",
				  "GGG" => "G"
				 );

sub translate
{
	my ($seq, $spos) = @_;
	$spos = 0 if !defined $spos;

	my $pep;
	for ( my $i = $spos; $i < length($seq); $i+=3 ) {
		last if $i+3 > length($seq); ## ignore non-triplet
		
		my $tri = substr $seq, $i, 3;
		$pep .= defined $CodonTable{$tri} ? $CodonTable{$tri} : "X";
	}
	return $pep;
}


1;
