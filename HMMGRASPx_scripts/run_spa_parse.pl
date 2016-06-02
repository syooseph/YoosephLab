# assume you are in the directory of saliva results i.e. $MY_SPA_DIR/Works/GRASPxp/Results/saliva_Butyrolactones

cat HMMGRASPx_out/recruited.list | awk '{print $2}' >graspxp.list 
perl ~/ScriptRepo/GetFastaSeqByID.pl --fasta=../../Data/saliva.faa --id=graspxp.list --out=graspxp.fa 
~/Tools/sfaspa-0.2.0/script/run_spa.pl -i graspxp.fa -s 1 -o SPA_graspxp/
~/Tools/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmsearch model.hmm SPA_graspxp/spa.fasta >graspxp.spa.results 
perl ../../Scripts/parse_spa_results.pl graspxp.fa graspxp.spa.results SPA_graspxp/
