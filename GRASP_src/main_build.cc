#include "extend_and_assemble.h"
#include "index_sample.h"
#include "gsa.h"
#include "sequence.h"
#include "guided_assemble.h"

#include <boost/program_options.hpp>
#include <iostream>
#include <string>

using namespace std;

typedef int AlignScoreType;
// define working directory and files to parse
static std::string working_directory = "WorkSpace";
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
  std::cout << "Usage: Build short_peptide_file" << std::endl;
  std::cout << "Use \'--help\' for more options" << std::endl;
  return;
}

int main(int argc, char** argv)	{
  InitDefaultParameters();
  
  boost::program_options::options_description desc("List of allowed options");
  desc.add_options()
      ("help", "print help message")
      ("short_peptides", boost::program_options::value<std::string>(&reads_mfasta_file), "multi_fasta file containing the short peptides")
      ("work_space", boost::program_options::value<std::string>(&working_directory)->default_value("WorkSpace"), "directory containing the indexing files")
      ("seed_len", boost::program_options::value<unsigned int>(&input_parameters.seed_len)->default_value(6), "length of the seeds")
      ("alphabet", boost::program_options::value<int>((int*) &input_parameters.alphabet)->default_value(4), "reduced alphabet to be used\n  0:ALL20\n  1:DSSP5\n  2:DSSP10\n  3:GBMR4\n  4:GBMR10\n  5:HSDM5\n  6:SDM6\n  7:MURPHY5\n  8:MURPHY10\n  9:TD5\n  10:TD10")
  ;
  boost::program_options::positional_options_description pos_opt;
  pos_opt.add("short_peptides", 1);
  
  boost::program_options::variables_map vm;
  boost::program_options::store(
      boost::program_options::command_line_parser(argc, argv).
      options(desc).positional(pos_opt).run(), vm
  );
  boost::program_options::notify(vm);
  
  if(vm.count("help"))  {
    std::cout << "Usage: Build short_peptide_file\n" << std::endl;
    std::cout << desc << std::endl; 
    return 0;
  }
  
  if(!vm.count("short_peptides"))  {
    std::cout << "The short peptides set which indexing is to be built on has not been specified." << std::endl;
    PrintUsage();
    return 0;
  }
  
  vector<string> files_in;
  files_in.resize(1);
  files_in[0] = reads_mfasta_file;
  int sample_num_reads = (unsigned int) seq::totalSequenceCount(files_in);
  //sample_tags_ = new char *[sample_num_reads_];
  char **sample_seqs = new char *[sample_num_reads];
  //cerr << "num of reads:  " << sample_num_reads_ << endl;
  char **tag_foo = NULL;
  seq::loadSequences(files_in, tag_foo, sample_seqs, SEQONLY);
  delete tag_foo;
  char **reversed_sample_seqs = new char *[sample_num_reads];
  //cerr << "num of reads:  " << sample_num_reads_ << endl;
  seq::reverseSequences(sample_num_reads, sample_seqs, reversed_sample_seqs);
  
  //cout << "Finish loading sequence" << endl;
  
  
  //if(working_directory.empty())  {
  //  working_directory = "./";
  //}
  //if(*(-- working_directory.end()) != '/')  {
  //  working_directory += "/";
  //}
  if(!boost::filesystem::exists(working_directory))  {
    boost::filesystem::path dir(working_directory);
    if(!boost::filesystem::create_directory(dir)) {
      cout << "GuidedAssemble::PrepareWorkSpace : Create working directory fails." << endl;
      exit(0);
    }  
  }
  GuidedAssemble func_caller;
  string read_stem = func_caller.GetFileStem(reads_mfasta_file);
  string sfa_file_name = working_directory + '/' + read_stem + ".sfa";
  string doc_file_name = working_directory + '/' + read_stem + ".doc";
  string lcp_file_name = working_directory + '/' + read_stem + ".lcp";
  string mcp_file_name = working_directory + '/' + read_stem + ".mcp";
  string gsa_file_name = working_directory + '/' + read_stem + ".gsa";
  
  //cout << "Finish infering file names" << endl;
  
  GSA* suffix_array = new GSA(sample_seqs, sample_num_reads, true);
  //cout << "Finish building suffix array" << endl;
  suffix_array->dump(
      sfa_file_name.c_str(), doc_file_name.c_str(), 
      lcp_file_name.c_str(), mcp_file_name.c_str(), 
      gsa_file_name.c_str()
  );
  //cout << "Finish dumping suffix array" << endl;
  
  sfa_file_name = working_directory + '/' + read_stem + ".re.sfa";
  doc_file_name = working_directory + '/' + read_stem + ".re.doc";
  lcp_file_name = working_directory + '/' + read_stem + ".re.lcp";
  mcp_file_name = working_directory + '/' + read_stem + ".re.mcp";
  gsa_file_name = working_directory + '/' + read_stem + ".re.gsa";
  GSA* reverse_suffix_array = new GSA(reversed_sample_seqs, sample_num_reads, true);
  //cout << "Finish building reverse suffix array" << endl;
  reverse_suffix_array->dump(
      sfa_file_name.c_str(), doc_file_name.c_str(), 
      lcp_file_name.c_str(), mcp_file_name.c_str(), 
      gsa_file_name.c_str()
  );
  //cout << "Finish dumping reverse suffix array" << endl;
  
  for(int i = 0; i < sample_num_reads; ++ i) {
    delete [] sample_seqs[i];
    delete [] reversed_sample_seqs[i];
  }
  delete [] sample_seqs;
  delete [] reversed_sample_seqs;
  
  string index_file_name = working_directory + '/' + read_stem + ".idx";
  string position_file_name = working_directory + '/' + read_stem + ".pos";
  
  IndexSample index_reads(input_parameters.seed_len, input_parameters.alphabet);
  if(!index_reads.IsCompatibleWithSetting(index_file_name, position_file_name))  {
    index_reads.WriteIndexFile(reads_mfasta_file, index_file_name, position_file_name);
  }
  
  return 0;
}
