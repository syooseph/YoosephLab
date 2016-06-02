#include "../include/kmer_unitcoder.h"

using namespace std;

KmerUnitcoder::KmerUnitcoder(const BioAlphabet &a, const int l) {
  alphabet_ = a;
  mer_len_ = l;
  // compute required bits
  for(num_bit_ = 1; num_bit_ <= 8 * sizeof(KmerUnitType); ++ num_bit_)  {
    if(pow(2, num_bit_) >= alphabet_.alphabet_size_)  {
      break;
    }
  }
  // check for valid setting
  CkValid();
  // bit-wise and with mask_ sets uncovered bits to 0
  mask_ = (KmerUnitType) 0;
  mask_ = ~mask_;
  for (int i = 0; i < mer_len_ * num_bit_; ++ i) {
    mask_ = mask_ << 1;
  }
  mask_ = ~mask_;
  // compute the signature for this BioAlphabet and unit k-mer length setting
  signature_ = (KmerUnitType) alphabet_.seq_type_;
  signature_ = signature_ << sizeof(KmerUnitType) * 4;
  KmerUnitType foo = (KmerUnitType) mer_len_;
  signature_ = signature_ | foo; 
  
  return;
}

KmerUnitcoder::~KmerUnitcoder(){
  return;
};

void KmerUnitcoder::CkValid()  {  
  if(num_bit_ * mer_len_ > 8 * sizeof(KmerUnitType))  {
    cout << "Error: cannot encode unit " << mer_len_ << "-mers within " << 8 * sizeof(KmerUnitType) << "-bit space, results will be corrupted" << endl;
    exit(1);
  }
  return;
}

KmerUnitType KmerUnitcoder::Encode(const char *s) {
  KmerUnitType c = 0;
  for(int i = 0; i < mer_len_; ++ i)  {
    KmerUnitType d = alphabet_.char_map_[s[i]];
    // bit-wise OR
    c = c << num_bit_;
    c = c | d;
  }
  return c;
}

int KmerUnitcoder::EncodeInt(const char *s) {
  int k = 0, n = alphabet_.GetSize();
  // check potential overflow
  if(pow(n, mer_len_) > 2147483646) {
    cout << "KmerUnicoder::Error: k-mer size too long, overflow detected. Abort." << endl;
    exit(1);
  }
  for(int i = 0; i < mer_len_; ++ i)  {
    int d = alphabet_.char_map_[s[i]];
    k = k * n + d;
  }
  return k;
}

int KmerUnitcoder::EncodeIntRight(const int num, const char c) {
  int n = alphabet_.GetSize();
  int d = alphabet_.char_map_[c];
  return (num % (int) pow(n, mer_len_ - 1)) * n + d;
}

std::string KmerUnitcoder::Decode(KmerUnitType ek)  {
  std::string s;
  for(int i = 0; i < mer_len_; ++ i)  {
    KmerUnitType d = pow(2, num_bit_) - 1;
    // bit-wise OR
    d = ek & d;
    ek = ek >> num_bit_;
    s += alphabet_.inv_char_map_[(int) d];
  }
  s = string(s.rbegin(), s.rend());
  return s;
}

KmerUnitType KmerUnitcoder::RightExt(const KmerUnitType ek, const char c) {
  KmerUnitType nek = ek;
  nek = nek << num_bit_;
  nek = nek | (KmerUnitType) alphabet_.char_map_[c];
  return nek & mask_;
}

KmerUnitType KmerUnitcoder::LeftExt(const KmerUnitType ek, const char c) {
  KmerUnitType nek = ek;
  nek = nek >> num_bit_;
  KmerUnitType d = (KmerUnitType) alphabet_.char_map_[c];
  d = d << num_bit_ * (mer_len_ - 1);
  nek = nek | d;
  return nek & mask_;
}


char KmerUnitcoder::HeadChar(const KmerUnitType ek) {
  KmerUnitType d = pow(2, num_bit_) - 1;
  d = d & (ek >> num_bit_ * (mer_len_ - 1));
  return alphabet_.inv_char_map_[(int) d];
}

char KmerUnitcoder::TailChar(const KmerUnitType ek)  {
  KmerUnitType d = pow(2, num_bit_) - 1;
  d = ek & d;
  return alphabet_.inv_char_map_[(int) d];
}

void KmerUnitcoder::PrintBits(const KmerUnitType ek) {
  bitset<8 * sizeof(KmerUnitType)> bek(ek);
  cout << bek << endl;
  return;
}
