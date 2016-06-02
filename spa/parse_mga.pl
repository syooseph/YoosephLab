#!/usr/bin/env perl

BEGIN { push @INC, `echo -n \$SPA` }

use strict;
use warnings;
use Set::Scalar;
#use Getopt::Long qw(GetOptionsFromArray);
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use io;
use seq;
use codon;

my ( $mga_file, $seq_file, $faa_file, $ffn_file, $help );
GetOptions( "mgafile=s" => \$mga_file,
			"seqfile=s" => \$seq_file,
			"orffile=s" => \$faa_file,
			"dnafile=s" => \$ffn_file,
			"help|h|?"  => \$help );

if ( !defined $mga_file || !defined $seq_file || $help ) {
	usage();
}

$faa_file = $mga_file . ".faa" unless defined $faa_file;
$ffn_file = $mga_file . ".ffn" unless defined $ffn_file;
	
my %Mga;
loadMetaGene();
proceed();

sub proceed
{
	my $orf = io::openOutput($faa_file);
	my $ffn = io::openOutput($ffn_file);
	my $sfh = io::openInput($seq_file);
	my ( $id, $dna );
	while ( <$sfh> ) {
		chomp;
		if ( /^>([^\s]+)/ ) {
			if ( defined $id && defined $Mga{$id} ) {
				extract( $orf, $ffn, $id, $dna );
			}
			$id = $1; $dna = "";
		} else {
			$dna .= $_;
		}
	}
	if ( defined $id && defined $Mga{$id} ) {
		extract( $orf, $ffn, $id, $dna );
	}
	close($sfh);
}

#sub translate
sub extract
{
	my ( $orf, $ffn, $id, $dna ) = @_;
	foreach my $entry ( @{$Mga{$id}} ) {
		#my ( $spos, $epos, $strand, $frame, $ends ) = @{$Mga{$id}};
		my ( $spos, $epos, $strand, $frame, $ends ) = @$entry;
		my $len = $epos-$spos+1;
# 		if ( $len <= 0 ) {
# 			print STDERR "Length Error:$id\t$len\n";
# 			print STDERR join(" ", @$entry), "\n";
# 		}
		my $sub = substr $dna, $spos-1, $len;
# 		if ( ($spos-1)+$len > length($dna) ) {
# 			print STDERR "Length Error:$id\t$spos\t$epos\t$len\n";
# 			print STDERR join(" ", @$entry), "\n";
# 		}
		$sub = seq::rc($sub) if $strand eq '-';
		
		my $pep = codon::translate($sub, $frame);
		my $nid = join("_", $id, $spos, $epos, $frame, $strand);

		print $ffn join("\n", ">$nid", $sub), "\n";
		print $orf join("\n", ">$nid", $pep), "\n";
	}
}

sub loadMetaGene
{
	my $id;
	my $mfh = io::openInput($mga_file);
	while ( <$mfh>) {
		chomp;
		next if /^# gc/;
		next if /^# self/;
# 		if ( /^# / ) {
# 			$id = $';
# 		} 
		## updated 07/18/2012
		if ( /^# ([^\s]+)/ ) {
			$id = $1;
		} 
		if ( /^gene/ ) {
			my @col = split;
			push @{$Mga{$id}}, [ @col[1..5] ];
		}
	}
	close($mfh);
}


sub usage
{
	print "Usage:$0 -m mga-file -s seq-file [-o orf-file] [-d ffn-file]\n";
	exit;
}

