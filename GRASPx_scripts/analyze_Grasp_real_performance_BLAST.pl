#!/usr/bin/perl -w
use strict;

my $sequence_dir = "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace";
my $script_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Scripts";
my $result_dir = "/usr/local/depot/projects/SPA/czhong/Grasp/Results";

chdir $result_dir or die "Cannot change directory: $!\n";
foreach(<LBPG*.fa.0.aln>)	{
	print "$_\n";
	my $aln_file = $_;
	$aln_file =~ /(.*)\.fa\.0\.aln/;
	my $l_tag = $1;
	# extracting the assembled sequence from the alignment file
	system "perl $script_dir/ConvertCertificateToMFASTA.pl $aln_file >$l_tag.asm.faa";
	# running blast between the query sequence and each of the assembled sequence
	system "perl $script_dir/run_reBLAST.pl $sequence_dir/$l_tag.fa $l_tag.asm.faa $l_tag.reBLAST.bst";
	# summarizing the results
	my $i;
	for($i = -10; $i <= 1; ++ $i)	{
		system "perl $script_dir/real_data_specificity.pl $l_tag.asm.faa $l_tag.reBLAST.bst 1e$i";
	}
}
