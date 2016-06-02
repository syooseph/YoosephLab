#include "sequence_build.h"
#include "reduced_alphabet.h"
#include "database_index.h"
#include "reachable_reads.h"
#include "scoring_function.h"
#include "read_alignment.h"
#include "greedy_assembly.h"
#include "assemble_extend.h"
#include "contig_refinement.h"
#include "remap.h"

#include "timer.h"

using namespace std;

int main(int argc, char** argv)  {
  
  double e_value = 10;
  int band_size = 20;
  double start_time = mytime();
  string file_in = "mix3fams.fa";
  SequenceBuild db_seq(file_in);
  db_seq.BuildSFADefault();
  db_seq.BuildKeyArrayDefault();
  SequenceBuild db_seq_rev(file_in);
  db_seq_rev.InPlaceReverse();
  db_seq_rev.BuildSFADefault();
  db_seq_rev.BuildKeyArrayDefault();
  DatabaseIndex db_dump(4, 6, 10);
  db_dump.BuildSeedmerMap(db_seq);
  unordered_map<std::string, std::list<std::string> > reduc_alph_map;
  db_dump.CreateReducedMap(reduc_alph_map);
  string out_file = file_in + ".rdm";
  db_dump.DumpReducedMap(reduc_alph_map, out_file);
  //reduc_alph_map.clear();
  //db_dump.LoadReducedMap(out_file, reduc_alph_map);
  unordered_map<string, list<PositionType> > fw_ext, re_ext;
  db_dump.CreateSeedExt(db_seq, db_seq_rev, fw_ext, re_ext);
  string out_file_seed_ext = file_in + ".sxt";
  db_dump.DumpSeedExt(fw_ext, re_ext, out_file_seed_ext);
  //db_dump.LoadSeedExt(out_file_seed_ext, fw_ext_seed, re_ext_seed);
  unordered_map<RIDType, list<OverlapType> > fw_read_ext, re_read_ext;
  db_dump.CreateReadExt(db_seq, db_seq_rev, fw_read_ext, re_read_ext);
  string out_file_read_ext = file_in + ".rxt";
  db_dump.DumpReadExt(fw_read_ext, re_read_ext, out_file_read_ext);
  //fw_read_ext.clear();
  //re_read_ext.clear();
  //db_dump.LoadReadExt(out_file_read_ext, fw_read_ext, re_read_ext);
  //return 0;
  
  /*****
  for(auto it = fw_read_ext.begin(); it != fw_read_ext.end(); ++ it) {
    cout << "INDEXING: READ EXTENSION GROUP:  " << it->first << endl;
    cout << "forward: ";
    for(auto itp = fw_read_ext[it->first].begin(); itp != fw_read_ext[it->first].end(); ++ itp) {
      cout << itp->rid << "," << (int) itp->len << "-";
    }
    cout << endl;
    cout << "reverse: ";
    for(auto itp = re_read_ext[it->first].begin(); itp != re_read_ext[it->first].end(); ++ itp) {
      cout << itp->rid << "," << (int) itp->len << "-";
    }
    cout << endl;
  }
  *****/
  double load_end_time = mytime();
  printElapsed(start_time, load_end_time, "load index");
  start_time = mytime();
  //return 0;
   
  ReachableReads candidate_selector(
      out_file, out_file_seed_ext, out_file_read_ext
  );
  
  ScoringFunction<int> score_scheme(PROTEIN, BLOSUM62, -1, -11);
  string query_seq = "MRAFIAIDVNESVRDSLVRAQDYIGSKEAKIKFVERENLHITLKFLGEITEEQAEEIKNILKKIAEKYKKHEVKVKGIGVFPNPNYIRVIWAGIENDEIIREMAREIEDELAKLGFKKEGNFVAHITLGRVKFVKDKLGLTMKLKELANEDFGSFVVDAIELKKSTLTPKGPIYETLARFELSE";
  //string query_seq = "MKYGIVLFPSKKLQDLANSYRKRYDPSYSLIPPHLTLRASFECAEEKADQLVSHLRNIAKESHPLVLKMTKYSSFAPVNNVIYIKAEPTEELKTLNEKLYTGVLAGEQEYNFVPHVTVGQNLSDDEHSDVLGQLKMQEVSHEEIVDRFHLLYQLENGSWTVYETFLLGRGE";
  //string query_seq = "MFNTTFANRTLPDLVEIQRASFCWFLNEGLAEEIQSFSPIVNYTGNLELHLFGDQYTLRYPKHNINECKRRDTTYSVQIYVPAQLINRETGVIKEQEVFIGDLPLMTDRGTFIINGAERVIVNQIVRSPGIYYKSELDKQGRRTYSSSLISNRGAWVKFETDRNDLVWVRIDKTRKIPAHVFLKAMGLSDTDIYNGLRHPEYLKKTFRFEGNYTTETALIQMYNKLRPGEPATVTGGQQLLYSRFFDPKRYDLGKVGRYKLNKKLNLSVPENVRVLTPQDTLAAIDYLINLKFEIGETDDIDHLGNRRVRSVGELLQNQVRIGLNRLERIIRERMTICDITSLTPNTLVNPKPIIASIREFFGSSQLSQFMDQTNPLAELTHKRRISALGPGGLNRDRAGFGVRDIHPSHYGRICPIETPEGPNAGLIGVLATHARINTYGFIEAPFFKVQDGQVYNHSQPIYLTADQEDKYRIAPGDITLDETNRIATKIVPIKYRQEFTTTKPNQVDFIAVSPIQVISIATSLIPFLEHDDANRALMGSNMQRQAVPLLYPESPLVGTGLEAQAARDSGMVVVSIEDGQVTFVSGDKICVTNKKGDEIAYYLQKYQRSNQDTCINQRPTVWLGEDVIEGQVIADGAATEGGELALGQNILVAYLPWEGYNYEDAFLINERLVYNDVYTSVHIEKYEIEARQTKLGSEEITRELPNVGEAALRKLDENGIIVIGSWVEAGDILIGKVTPKGESDQPPEGKLLRAIFGEKARDVRDTSLRVPNGGRGRILDVRIFTREKGDELPTGANIVIRVYIAQSRKIQVGDKMAGRHGNKGIISRILPRQDMPYLPDGTPVDLVLNPLGVPSRMNVGQIFECLLGLAAENLNKRFKITPFDEMHGAEASRVLVNEKLNEAKIKTGENWLFDLRHPGKITLYDGRTGEAFDNPVTIGVSYMLKLVHLVDDKIHARSTGPYSLVTQQPLGGRAQHGGQRLGEMEVWALEAFGASYTLQELLTVKSDDMQGRNETLNAIVKGKPIPRPGTPESFKVLMRELQSLGLDIGAYKIENLPDGQTRGIEVDLMMNYQQSRLFKPLYESMQTKNNENLFL";
  map<int, list<SeedType> > candidate_seeds;
  candidate_selector.SelectSeeds(score_scheme, query_seq, candidate_seeds);
  
  double select_end_time = mytime();
  printElapsed(start_time, select_end_time, "select seeds");
  start_time = mytime();
  
  /******
  for(auto it = candidate_seeds.rbegin(); it != candidate_seeds.rend(); ++ it) {
    cout << "alignment score: " << it->first << endl;
    for(auto it_r = it->second.begin(); it_r != it->second.end(); ++ it_r) {
      cout << it_r->seed_seq << "." << it_r->q_pos << ",";
    }
    cout << endl;
  }
  ******/
  
  //map<int, list<ReadType> > seed_reads;
  //candidate_selector.GetSeedReads(db_seq, candidate_seeds, seed_reads);
  
  /******
  for(auto it = seed_reads.rbegin(); it != seed_reads.rend(); ++ it) {
    cout << "****** score: " << it->first << endl;
    for(auto it_r = it->second.begin(); it_r != it->second.end(); ++ it_r) {
      cout << it_r->rid << "  " << it_r->r_pos << " " << db_seq.GetSequence(it_r->rid).substr((int) it_r->r_pos, 6) << endl;
    }
  }
  ******/
  
  //list<RIDType> candidate_reads;
  //candidate_selector.CollectNeighbourReads(db_seq, 10, candidate_seeds, candidate_reads);
  /**
  cout << "Recruited reads: " << endl;
  for(auto it = candidate_reads.begin(); it != candidate_reads.end(); ++ it) {
    cout << it->rid << "  " << db_seq.GetHeader(it->rid) << endl;
  }
  **/
  //double collect_end_time = mytime();
  //printElapsed(start_time, collect_end_time, "collect reads");
  //start_time = mytime();
  
  //list<ContigType> contigs;
  //AssembleExtend guided_assembly;
  //guided_assembly.AssembleAllContigs(
  //    query_seq, db_seq, candidate_selector, score_scheme, band_size, 
  //    seed_reads, candidate_reads, e_value * 100, contigs
  //);
  //int n = 0;
  //for(auto it = contigs.begin(); it != contigs.end(); ++ it) {
    //cout << "Final contig:  " << it->bit_score << "(" << it->score << ")" << " " << it->sequence  << endl;
    //cout << ">contig_old_" << n << ":" << it->q_begin << "-" << it->q_end << endl << it->sequence << endl;
  //  ++ n;
  //}
  
  //double extend_end_time = mytime();
  //printElapsed(start_time, extend_end_time, "extend contigs");
  //start_time = mytime();
  
  //list<ContigType> refined_contigs;
  //ContigRefinement recalibrate_contigs;
  //recalibrate_contigs.RefineContigs(
  //    query_seq, db_seq, score_scheme, band_size, e_value, 10, contigs, refined_contigs
  //);
  //n = 0;
  //for(auto it = refined_contigs.begin(); it != refined_contigs.end(); ++ it) {
    //cout << "Final contig:  " << it->bit_score << "(" << it->score << ")" << " " << it->sequence  << endl;
    //cout << ">contig_new_" << n << ":" << it->q_begin << "-" << it->q_end << endl << it->sequence << endl;
  //  ++ n;
  //}
  
  //double refine_end_time = mytime();
  //printElapsed(start_time, refine_end_time, "refine contigs");
  //start_time = mytime();
  
  //**********
  //ReMap re_aln;
  //list<MapReadType> identified_reads;
  //re_aln.MapToContig(db_seq, refined_contigs, 30, 0.9, 1, false, identified_reads);
  //for(auto it = identified_reads.begin(); it != identified_reads.end(); ++ it) {
  //  cout << it->rid << endl;
  //}
  //**********/
  
  //double remap_end_time = mytime();
  //printElapsed(start_time, remap_end_time, "re-align reads");
  //start_time = mytime();
  
  //ReadAlignment map_reads;
  //map_reads.ComputeAlignmentScore(query_seq, db_seq, score_scheme, 10, candidate_reads);
  //map_reads.SortReadsOnEvalue(candidate_reads);
  
  
  
  /*
  GreedyAssembly ref_assembler;
  unordered_map<RIDType, list<ReadListIterType> > rid_link;
  ref_assembler.BuildRIDLink(candidate_reads, rid_link);
  int n = 0;
  for(auto it = candidate_reads.begin(); it != candidate_reads.end(); ++ it) {
    std::deque<ReadConnectType> contig_path;
    ref_assembler.GreedyExtend(candidate_selector, rid_link, *it, 1, contig_path);
    //cout << "********** seed: " << it->rid << endl;
    //for(auto it_r = contig_path.begin(); it_r != contig_path.end(); ++ it_r) {
    //  cout << it_r->rid << ","; 
    //}
    //cout << endl;
    string spell_seq = ref_assembler.SpellContigSequence(db_seq, contig_path);
    cout << ">contig_" << n << endl << spell_seq << endl;
    ++ n;
  }
  */
  
  //vector<deque<ReadConnectType> > pre_contigs;
  //unordered_map<int, list<int> > fw_connect;
  //unordered_map<int, list<int> > re_connect;
  //candidate_selector.ConnectReads(candidate_reads, pre_contigs, fw_connect, re_connect);
  
  //for(auto it = pre_contigs.begin(); it != pre_contigs.end(); ++ it) {
  //  cout << "********** contig size:  " << it->size() << endl; 
  //  for(auto it_r = it->begin(); it_r != it->end(); ++ it_r) {
  //    cout << it_r->rid << ",";
  //  }
  //  cout << endl;
  //}
  
  /*****
  for(auto it_fw = fw_ext.begin(); it_fw != fw_ext.end(); ++ it_fw) {
    cout << "********** " << it_fw->first << endl;
    for(auto it_rf = fw_ext[it_fw->first].begin(); it_rf != fw_ext[it_fw->first].end(); ++ it_rf) {
      cout << (unsigned int) it_rf->rid << "," << (unsigned int) it_rf->pos << ";";
    }
    cout << endl;
    for(auto it_rf = re_ext[it_fw->first].begin(); it_rf != re_ext[it_fw->first].end(); ++ it_rf) {
      cout << (unsigned int) it_rf->rid << "," << (unsigned int) it_rf->pos << ";";
    }
    cout << endl;
  } 
  *****/
  /*****
  for(auto it = reduc_alph_map.begin(); it != reduc_alph_map.end(); ++ it) {
    cout << it->first << endl;
    for(auto it_p = it->second.begin(); it_p != it->second.end(); ++ it_p) {
      cout << *it_p << ";";
    }
    cout << endl;
  }
  *****/
  
  return 0;
}

