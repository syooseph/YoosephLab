#include <inttypes.h>
#include <math.h>
#include <unistd.h>
#include <assert.h>

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <tuple>
#include <string>
#include <list>

// The class is used to encode all k-mer positions of the sample,
// and put them into the harddist.
// The inputs required are the k-mer size and the alphabet to be used.
// Sample usage:
//  IndexSample seq_encoder(7, GBMR10); // 7 is the k-mer size,
//                                      // GBMR10 is a supported alphabet
//  // The following command writes all kmer info into harddisk
//  seq.encoder.WriteIndexFile(
//              "my_seq.fasta",         // my_seq.fasta is a multi-fasta file 
//              "my_seq.fasta.idx",     // my_seq.fasta.idx is the output index file
//              "my_seq.fasta.pos");    // my_seq.fasta.pos is the output position file
//  // The following command read all kmer position in query_seq into position_holder
//  std::unordered_map<KmerType, list<PositionType> > position_holder;
//  seq.ReadKmerIndex(my_seq.fasta.idx,
//                    my_seq.fasta.pos,
//                    query_seq,
//                    position_holder);

#ifndef _INDEX_SAMPLE_H_
#define _INDEX_SAMPLE_H_

typedef int64_t KmerType;
typedef unsigned int RIDType;
typedef unsigned char POSType;

struct PositionType {
  // the 32bit read ID
  RIDType rid;
  // the 8bit read position
  POSType pos;
  
  PositionType& operator= (const PositionType& source)  {
    this->rid = source.rid;
    this->pos = source.pos;
    return *this;
  }
};

inline std::string ReverseShortSequence(const std::string seq)  {
  std::string rev_seq = std::string(seq.rbegin(), seq.rend());
  return rev_seq;
}

// a list of supported alphabets
enum Alphabet {ALL20, DSSP5, DSSP10, GBMR4, GBMR10, HSDM5, SDM6, MURPHY5, MURPHY10, TD5, TD10};

class IndexSample  {
  
 public:
 
  IndexSample(void);
  IndexSample(unsigned int in_kmer_size, enum Alphabet in_alphabet);
  ~IndexSample(void);
  void WriteIndexFile(std::string in_file, std::string out_index_file, std::string out_position_file);
  void ReadKmerIndex(
    std::string in_index_file, 
    std::string in_position_file, 
    std::string query_seq, 
    std::unordered_map<KmerType, std::list<PositionType> >& kmer_locations
  );
  void ReadKmerIndex(
    std::string in_index_file, 
    std::string in_position_file, 
    std::vector<std::string> query_seq, 
    std::unordered_map<KmerType, std::list<PositionType> >& kmer_locations
  );
  void BuiltParamLoader(
    const std::string& in_index_file, const std::string& in_position_file,
    enum Alphabet& alp_in_use, unsigned int& ksize_in_use
  );
  int AccessAlphabetMap(char c);
  // the conversion functions
  inline PositionType encode_read_position(const RIDType read_id, const POSType position) {
    // ensure the input is not out of range
    PositionType encoded;
    assert(sizeof(encoded.rid) == 4 && sizeof(encoded.pos) == 1);
    assert(log2(read_id) < 32);
    assert(log2(position) < 8);
    // the rid field corresponds to the read_id field
    encoded.rid = read_id;
    encoded.pos = position;
    return encoded;
  }
  inline void decode_read_position(const PositionType encoded_position, 
      unsigned int &read_id, unsigned int &position)  
  {
    assert(sizeof(encoded_position.rid) == 4 && sizeof(encoded_position.pos) == 1);
    read_id = position = 0;
    read_id = (unsigned int) encoded_position.rid;
    position = (unsigned int) encoded_position.pos;
    return;
  }
  inline KmerType encode_kmer(const std::string kmer_seq)  {
    if(KmerSize != kmer_seq.length())  {
      std::cout << KmerSize << std::endl;
      std::cout << kmer_seq << std::endl;
    }
    assert(KmerSize == kmer_seq.length());
    KmerType encoded_kmer = 0;
    unsigned int i, c;
    for(i = 0; i < KmerSize; ++ i)  {
      c = AlphabetMap[kmer_seq[i]];
      encoded_kmer |= c;
      encoded_kmer <<= EncodeCharacterBits;
    }
    return encoded_kmer;
  }
  
  bool IsCompatibleWithSetting(const std::string& idx_file, const std::string& pos_file);
  
  // Obsolute function ReconsolidateIndex, should not be used
  // void ReconsolidateIndex(
  //  std::string in_index_file, std::string in_position_file, 
  //  std::string out_index_file, std::string out_position_file 
  // );

 protected:
  int AlphabetSize;
  unsigned int KmerSize;
  int EncodeCharacterBits;
  enum Alphabet CurrentAlphabet;
  
  unsigned int KmerHashEntries;
  unsigned int KmerHashMaxLength;
  
  std::map<char, int> AlphabetMap;
  std::unordered_map<KmerType, std::list<PositionType> > KmerHash;
  
  void init_alphabet(enum Alphabet a);
  
  void assign_alphabet_ALL20();
  void assign_alphabet_DSSP5();
  void assign_alphabet_DSSP10();
  void assign_alphabet_GBMR4();
  void assign_alphabet_GBMR10();
  void assign_alphabet_HSDM5();
  void assign_alphabet_SDM6();
  void assign_alphabet_MURPHY5();
  void assign_alphabet_MURPHY10();
  void assign_alphabet_TD5();
  void assign_alphabet_TD10();
  
  double get_log2_sys_mem();

  void check_setup(void);
  inline void dump_to_HD(const unsigned int cutoff, std::ostream& idx_fh, std::ostream& pos_fh);
  inline void update_position(unsigned int read_id, std::string seq, std::ostream& idx_fh, std::ostream& pos_fh);
  void write_format_identifier(std::string identifier, std::ostream& file);
  std::string read_format_identifier(std::istream& file, enum Alphabet& current_alphabet, unsigned int& kmer_size);
  std::string read_format_identifier(std::istream& file);
  
  
  void get_position_file_offsets(
    enum Alphabet alphabet,
    unsigned int kmer_size,
    std::string query_seq, 
    std::istream& idx_fh, 
    std::map<std::streampos, std::pair<KmerType, std::streampos> >& file_offsets
  );
  
  void get_position_file_offsets(
    enum Alphabet alphabet,
    unsigned int kmer_size,
    std::vector<std::string> query_seq, 
    std::istream& idx_fh, 
    std::map<std::streampos, std::pair<KmerType, std::streampos> >& file_offsets
  );
  
  void get_kmer_positions(
    std::istream& pos_fh,
    // the input start and end trunks to be read
    const std::map<std::streampos, std::pair<KmerType, std::streampos> >& file_offsets,
    // the output tells all positions for each k-mer
    std::unordered_map<KmerType, std::list<PositionType> >& kmer_positions
  );
  
  void load_index_file(
    std::string in_index_file, 
    std::unordered_map<KmerType, std::list< std::pair<std::streampos, std::streampos> > >& kmer_positions
  );
  
  inline void load_position_per_kmer(
    std::istream &pos_fh, 
    std::unordered_map<KmerType, std::list< std::pair<std::streampos, std::streampos> > >::const_iterator it,
    std::list<PositionType>& loaded_positions
  );
  
  bool ContainSpecialChar(std::string mer);
};

#endif
