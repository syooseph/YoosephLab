##==============
## Generate kmer
##============== 
  ./kmer
  e.g. ./kmer -k 6 ../data/reca.faa -o /tmp/reca.k6.list

##======================
## k-mer to read mapping
##======================
1. write a map in no order
   ./kr
   e.g. nohup time ./kr -n /usr/local/scratch/yyang/fgs/fgs.read.names.v2 -d 1 -t /usr/local/scratch/yyang/fgs/kr1 /usr/local/scratch/yyang/*faa >& $SCRATCH/fgs/kr1.log &
   - read-name to read-id mapping
   - 24 read-id files (1 file for 1 AA), which splitted from 1st file
2. combine maps
   ./combine.sh
   e.g. nohup time ./combine.sh /usr/local/scratch/yyang/fgs/kr1 6 >& /usr/local/scratch/yyang/fgs/combine.log &

##==============
## Run assembler
##==============
  ./spa
  e.g. ./spa -k 6 -m 3 -r 10 --no-transform --dump -t ../data/reca/kr -i ../data/reca.6mer.list >& /tmp/zzz7
