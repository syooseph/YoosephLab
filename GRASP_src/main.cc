#include "extend_and_assemble.h"
#include "index_sample.h"
#include "gsa.h"
#include "sequence.h"
#include "guided_assemble.h"

#include <boost/program_options.hpp>
#include <iostream>
#include <string>

typedef int AlignScoreType;
// define working directory and files to parse
static std::string working_directory = "WorkSpace";
static std::string query_file;
static std::string reads_mfasta_file;
// setting default parameters
static AssembleParameterType<AlignScoreType> input_parameters;

void InitDefaultParameters()  {
  input_parameters.num_threads = 1;
  input_parameters.seed_len = 6;
  input_parameters.alphabet = GBMR10;
  input_parameters.gap_open = -11;
  input_parameters.gap_extension = -1;
  input_parameters.scoring_matrix = BLOSUM62;
  input_parameters.band_size = 40;
  input_parameters.evalue_cutoff = 10;
  input_parameters.n_back_check = 10;
  return;
}

void PrintUsage() {
  std::cout << "Usage: grasp query_file short_peptide_file" << std::endl;
  std::cout << "Use \'--help\' for more options" << std::endl;
  return;
}

int main(int argc, char** argv)	{
  InitDefaultParameters();
  
  boost::program_options::options_description desc("List of allowed options");
  desc.add_options()
      ("help", "print help message")
      ("query", boost::program_options::value<std::string>(&query_file), "multi-fasta file containing the query sequences")
      ("short_peptides", boost::program_options::value<std::string>(&reads_mfasta_file), "multi_fasta file containing the short peptides")
      ("work_space", boost::program_options::value<std::string>(&working_directory), "directory containing the indexing files")
      ("seed_len", boost::program_options::value<unsigned int>(&input_parameters.seed_len)->default_value(6), "length of the seeds")
      ("alphabet", boost::program_options::value<int>((int*) &input_parameters.alphabet)->default_value(4), "reduced alphabet to be used\n  0:ALL20\n  1:DSSP5\n  2:DSSP10\n  3:GBMR4\n  4:GBMR10\n  5:HSDM5\n  6:SDM6\n  7:MURPHY5\n  8:MURPHY10\n  9:TD5\n  10:TD10")
      ("scoring_matrix", boost::program_options::value<int>((int*) &input_parameters.scoring_matrix)->default_value(0), "scoring matrix to be used\n  0:BLOSUM62\n  1:BLOSUM80\n  2:BLOSUM90\n  3:BLOSUM50\n  4:BLOSUM45\n  5:PAM250\n  6:PAM70\n  7:PAM30")
      ("gap_extension", boost::program_options::value<AlignScoreType>(&input_parameters.gap_extension)->default_value(-1), "gap extension penalty")
      ("gap_open", boost::program_options::value<AlignScoreType>(&input_parameters.gap_open)->default_value(-11), "gap open penalty")
      ("band_size", boost::program_options::value<unsigned int>(&input_parameters.band_size)->default_value(40), "size for banded sequence alignment")
      ("e_value", boost::program_options::value<double>(&input_parameters.evalue_cutoff)->default_value(10), "e-value cutoff")
      ("back_check", boost::program_options::value<unsigned int>(&input_parameters.n_back_check)->default_value(10), "minimum overlap with the existing path while being extended")
      ("num_threads", boost::program_options::value<unsigned int>(&input_parameters.num_threads)->default_value(1), "number of threads")
  ;
  boost::program_options::positional_options_description pos_opt;
  pos_opt.add("query", 1);
  pos_opt.add("short_peptides", 1);
  
  boost::program_options::variables_map vm;
  boost::program_options::store(
      boost::program_options::command_line_parser(argc, argv).
      options(desc).positional(pos_opt).run(), vm
  );
  boost::program_options::notify(vm);
  
  if(vm.count("help"))  {
    std::cout << "Usage: grasp query_file short_peptide_file\n" << std::endl;
    std::cout << desc << std::endl; 
    return 0;
  }
  
  if(!vm.count("query") || !vm.count("short_peptides"))  {
    std::cout << "The query sequence or the short peptides is not specified." << std::endl;
    PrintUsage();
    return 0;
  }
  
  GuidedAssemble align_and_assemble(
      working_directory, query_file, reads_mfasta_file, input_parameters
  );
  align_and_assemble.PrepareRunEnvironment();
  
  //return 0;  

  std::unordered_map<RIDType, ReadMapType> mapped_reads;
  align_and_assemble.AssembleAndAlign(mapped_reads);
  std::vector<std::string> read_tags;
  align_and_assemble.GetReadTags(read_tags);
  for(auto it = mapped_reads.begin(); it != mapped_reads.end(); ++ it) {
    //double bit_score = score_scheme.ComputeBitScore((it->second).assembled_score);
    std::cout << read_tags[static_cast<int>(it->first)] << "  " << it->second.best_evalue << std::endl;
    //std::cout << (unsigned int) it->first << "  " << bit_score << std::endl;
  }
  //std::cout << "num of reads mappings:  " << mapped_reads.size() << std::endl;
  
  return 0;
}
