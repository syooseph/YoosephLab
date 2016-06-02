#include "guided_assemble.h"
#include "seq_align.h"
#include <cstring>

typedef int AlignScoreType;

class RealignmentFunctor  {
 public:
  RealignmentFunctor(
      void *in_saa_obj, const std::string& in_query_seq, const std::string& in_target_seq, 
      AssembledPathType* in_path, std::list<CertificateType<AlignScoreType> > *in_certificate_holder,
      LockType* in_mutex_holder
  ) :  saa_obj_(in_saa_obj),
      query_seq_(in_query_seq),
      target_seq_(in_target_seq),
      path_(in_path),
      certificate_holder_(in_certificate_holder),
      mutex_holder_(in_mutex_holder)  
  {
    return;
  }
  ~RealignmentFunctor() {return;}
  void run();
 private:
  void *saa_obj_;
  std::string query_seq_;
  std::string target_seq_;
  AssembledPathType* path_;
  std::list<CertificateType<AlignScoreType> > *certificate_holder_;
  LockType* mutex_holder_;
};

void RealignmentFunctor::run()  {
  GuidedAssemble* func_caller = (GuidedAssemble*) saa_obj_;
  SeqAlign<AlignScoreType> re_alignment(
      target_seq_, query_seq_, func_caller->score_scheme_, SEMIGLOBAL
  );
  re_alignment.Align();
  re_alignment.TraceBack();
  std::unordered_map<int, int> nuc_matching;
  CertificateType<AlignScoreType> certificate;
  re_alignment.GetAlignment(
      nuc_matching, certificate.alignment_target, 
      certificate.alignment_query, certificate.alignment_symbol
  );
  certificate.alignment_score = re_alignment.GetBestScore();
  certificate.alignment_regions = re_alignment.GetAlignmentRegion();
  for(auto it = path_->mapped_locations.begin(); it != path_->mapped_locations.end(); ++ it) {
    if(it->second >= certificate.alignment_regions[0] && 
        it->second + 31 <= certificate.alignment_regions[1])  {
      certificate.mapped_locations.push_back({it->first, it->second - certificate.alignment_regions[2]});
    }
  }
  ULock w_lock(*mutex_holder_);
  certificate_holder_->push_back(certificate);
  w_lock.unlock();
  return;
}


