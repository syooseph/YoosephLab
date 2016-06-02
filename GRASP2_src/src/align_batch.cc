#include "../include/align_batch.h"

using namespace std;

int AlignBatch::Max3(int a, int b, int c) {
  int t = a > b ? a : b;
  return t > c ? t : c;
}

int AlignBatch::Min3(int a, int b, int c) {
  int t = a < b ? a : b;
  return t < c ? t : c;
}

bool SortFullAlnType(const FullAlnType &a, const FullAlnType &b)  {
  if((a.score > b.score) || 
    (a.score == b.score && a.alignment.first.length() < b.alignment.first.length())
  ) return true;
  return false;
}

void AlignBatch::MultiAlignGlobal(
    const int threads, char *sp, const int n, char **sgroup, 
    const int band, ScoringProt &f, std::vector<int> &r
) {
  int i, j, a, b, l = strlen(sp);
  // define the effective band: if the sequence in the target is too long
  // we cannot guarantee that the banding is enough for covering the entire target sequence
  // need to adjust it
  vector<int> eff_band(threads, 0);
  // setting the effective length through refering to the longest sequence int the set
  vector<int> range;
  range.push_back(-1);
  int chunk_size = n / threads;
  for(i = 0; i < threads - 1; ++ i) range.push_back(range.back() + chunk_size);
  range.push_back(n - 1);
  // finding the longest length of the chunk and recompute the effective banding
  for(i = 0; i < threads; ++ i) {
    for(j = range[i] + 1; j <= range[i + 1]; ++ j) {
      // now eff_band[i] stores the lagest difference between the sequence in the chunk
      // and the query sequence
      a = abs(l - (int) strlen(sgroup[j]));
      eff_band[i] = eff_band[i] > a ? eff_band[i] : a;
    }
    // increase the band size by "band"
    eff_band[i] *= 2;  eff_band[i] += band;
  } 
  // execute the alignment in multi-threaded mode
  vector<vector<int> > r_multi(threads);
  #pragma omp parallel num_threads(threads)
  {
    #pragma omp for
    for(i = 0; i < threads; ++ i) {
      AlignGlobal(
          sp, range[i + 1] - range[i], &sgroup[range[i] + 1], 
          eff_band[i], f, r_multi[i]
      );
    }
    #pragma omp taskwait
  }
  // merge the alignment score results
  r.resize(n); a = 0;
  for(i = 0; i < threads; ++ i) {
    for(j = 0; j < r_multi[i].size(); ++ j) r[a ++] = r_multi[i][j];
  }
  return;
}

void AlignBatch::MultiAlignGlobal(
    const int threads, std::string &sp, const int n, std::vector<std::string> &sgroup, 
    const int band, ScoringProt &f, std::vector<int> &r
) {
  int i, j, a, b, l = sp.length();
  // define the effective band: if the sequence in the target is too long
  // we cannot guarantee that the banding is enough for covering the entire target sequence
  // need to adjust it
  vector<int> eff_band(threads, 0);
  // setting the effective length through refering to the longest sequence int the set
  vector<int> range;
  range.push_back(-1);
  int chunk_size = n / threads;
  for(i = 0; i < threads - 1; ++ i) range.push_back(range.back() + chunk_size);
  range.push_back(n - 1);
  // finding the longest length of the chunk and recompute the effective banding
  for(i = 0; i < threads; ++ i) {
    for(j = range[i] + 1; j <= range[i + 1]; ++ j) {
      // now eff_band[i] stores the lagest difference between the sequence in the chunk
      // and the query sequence
      a = abs(l - (int) sgroup[j].length());
      eff_band[i] = eff_band[i] > a ? eff_band[i] : a;
    }
    // increase the band size by "band"
    eff_band[i] *= 2;  eff_band[i] += band;
  } 
  // execute the alignment in multi-threaded mode
  vector<vector<int> > r_multi(threads);
  #pragma omp parallel num_threads(threads)
  {
    #pragma omp for
    for(i = 0; i < threads; ++ i) {
      vector<string> schunk(sgroup.begin() + range[i] + 1, sgroup.begin() + range[i + 1] + 1);
      AlignGlobal(
          sp, range[i + 1] - range[i], schunk, 
          eff_band[i], f, r_multi[i]
      );
    }
    #pragma omp taskwait
  }
  // merge the alignment score results
  r.resize(n); a = 0;
  for(i = 0; i < threads; ++ i) {
    for(j = 0; j < r_multi[i].size(); ++ j) r[a ++] = r_multi[i][j];
  }
  return;
}

void AlignBatch::MultiAlignGlobalPairwise(
    const int threads, 
    std::vector<std::string> &sgroup1, std::vector<std::string> &sgroup2,
    const int band, ScoringProt &f, std::vector<int> &r
) {
  if(sgroup1.size() != sgroup2.size())  {
    cout << "Warning: AlignBatch::MultiAlignLocalPairwise: number of sequences in each group do not match!!!, abort." << endl;
    exit(1);
  }
  int i, j, n = sgroup1.size();
  // computing the effective band size
  vector<int> eff_band(n, band);
  for(i = 0; i < n; ++ i) {
    eff_band[i] = abs((int) sgroup1[i].size() - (int) sgroup2[i].size()) * 2 + band;
  }
  r.resize(n);
  // compute the alignments
  #pragma omp parallel num_threads(threads)
  {
    #pragma omp for
    for(i = 0; i < n; ++ i) {
      r[i] = AlignGlobalPairwise(sgroup1[i], sgroup2[i], band, f);
    }
    #pragma omp taskwait
  }
  return;
}

void AlignBatch::AlignGlobal(
    char *sp, const int n, char **sgroup, 
    const int band, ScoringProt &f, 
    std::vector<int> &r
) {
  int lsp = strlen(sp);
  if(lsp <= 0 || n <= 0 || band < 2)  return;
  int i, j, k, x, t, a, b, c, s, rl;
  int g = f.gap_open_, e = f.gap_ext_;
  // declaring the checking matrices
  bool *filled = new bool [n]; bool *finished = new bool [n]; r.resize(n, -MAX);
  // mm: matrix match; mi: matrix insert; md: matrix delete
  int **mm = new int* [n], **mi = new int* [n], **md = new int*[n];
  for(i = 0; i < n; ++ i) {
    filled[i] = finished[i] = false;
    mm[i] = new int [band]; mi[i] = new int [band]; md[i] = new int [band];
  }
  // initalizing the scores
  int bh = band / 2;
  for(i = 0; i < bh; ++ i) mm[0][i] = mi[0][i] = md[0][i] = -MAX;
  mm[0][bh] = 0; mi[0][bh] = md[0][bh] = g;
  for(i = bh + 1; i < band; ++ i)  {
    mm[0][i] = mi[0][i] = mi[0][i - 1] + e;
    md[0][i] = mm[0][i] + g;
  }
  // note that these two entries are never used; only set for easier computation of other entries
  int z = sizeof(int) * band;
  for(i = 1; i < n; ++ i) {
    memcpy(mm[i], mm[0], z); memcpy(mi[i], mi[0], z); memcpy(md[i], md[0], z); 
  }
  
  int *s1_leftover = new int [n], *s2_leftover = new int [n];
  for(i = 0; i < n; ++ i) {
    s1_leftover[i] = lsp; s2_leftover[i] = strlen(sgroup[i]) - (band - bh) + 1;
  } 
  // perform alignment
  for(i = 1; i <= lsp; ++ i) {   // for each iteration or query sequence position
    for(k = 0; k < n; ++ k) { // for each sequence in the target
      if(finished[k]) continue;
      rl = strlen(sgroup[k]);
      for(x = 0; x < band; ++ x) { // for each band region in the specific target sequence
        j = i - bh + x; // computing the original index of the position
        if(j < 0) continue;
        if(j > rl) {
          if(x == 0)  {finished[k] = true;} // no need to look at the sequence any more 
          else  {
            filled[k] = true; 
            r[k] = Max3(mm[k][x - 1], mi[k][x - 1], md[k][x - 1]);
          } // record score anyway (global alignment)
          break;  // no need to continue computing for this band any more
        }
        if(j >= 1) s = f.GetMatchScore(sp[i - 1], sgroup[k][j - 1]);
        // computing order: Insertion matrix (mi) -> Match matrix (mm) -> Deletion matrix (md)
        if(j == 0)  {
          mm[k][x] = md[k][x] = md[k][x + 1] + e; mi[k][x] = mm[k][x] + g;
        } else if(x == band - 1) {
          a = mi[k][x - 1] + e; b = mm[k][x - 1] + g + e;
          mi[k][x] = a > b ? a : b;
          c = mm[k][x] + s;
          mm[k][x] = a > c ? a : c; 
          md[k][x] = -MAX;
        } else  {
          if(x <= 0) {a = -MAX;} else {a = mi[k][x - 1] + e;}
          if(x <= 0) {b = -MAX;} else {b = mm[k][x - 1] + g + e;}
          mi[k][x] = a > b ? a : b;
          b = mm[k][x + 1] + g + e; c = md[k][x + 1] + e;
          md[k][x] = b > c ? b : c;
          t = mm[k][x] + s;
          mm[k][x] = mi[k][x] > md[k][x] ? mi[k][x] : md[k][x];
          mm[k][x] = t > mm[k][x] ? t : mm[k][x]; 
        }
      }
      if(!finished[k]) {
        -- s1_leftover[k]; -- s2_leftover[k];
      }
    }
  }
  // check for the last time for scores that have not been filled
  // if not filled, fill it with the score that has been latestly computed
  for(k = 0; k < n; ++ k) {
    if(!filled[k]) r[k] = Max3(mm[k][band - 1], mi[k][band - 1], md[k][band - 1]);
    if(s1_leftover[k] > 0) r[k] += g + s1_leftover[k] * e;
    if(s2_leftover[k] > 0) r[k] += g + s2_leftover[k] * e;
  }
  // collect memory
  for(i = 0; i < n; ++ i) {delete [] mi[i]; delete [] mm[i]; delete [] md[i];}
  delete [] mi; delete [] mm; delete [] md;
  delete [] filled; delete [] finished;
  delete [] s1_leftover; delete [] s2_leftover;
  return;
}

void AlignBatch::AlignGlobal(
    std::string &sp, const int n, std::vector<std::string> &sgroup, 
    const int band, ScoringProt &f, 
    std::vector<int> &r
) {
  int lsp = sp.length();
  if(lsp <= 0 || n <= 0 || band < 2)  return;
  int i, j, k, x, t, a, b, c, s, rl;
  int g = f.gap_open_, e = f.gap_ext_;
  // declaring the checking matrices
  bool *filled = new bool [n]; bool *finished = new bool [n]; r.resize(n, -MAX);
  // mm: matrix match; mi: matrix insert; md: matrix delete
  int **mm = new int* [n], **mi = new int* [n], **md = new int*[n];
  for(i = 0; i < n; ++ i) {
    filled[i] = finished[i] = false;
    mm[i] = new int [band]; mi[i] = new int [band]; md[i] = new int [band];
  }
  // initalizing the scores
  int bh = band / 2;
  for(i = 0; i < bh; ++ i) mm[0][i] = mi[0][i] = md[0][i] = -MAX;
  mm[0][bh] = 0; mi[0][bh] = md[0][bh] = g;
  for(i = bh + 1; i < band; ++ i)  {
    mm[0][i] = mi[0][i] = mi[0][i - 1] + e;
    md[0][i] = mm[0][i] + g;
  }
  // note that these two entries are never used; only set for easier computation of other entries
  int z = sizeof(int) * band;
  for(i = 1; i < n; ++ i) {
    memcpy(mm[i], mm[0], z); memcpy(mi[i], mi[0], z); memcpy(md[i], md[0], z); 
  }
  
  int *s1_leftover = new int [n], *s2_leftover = new int [n];
  for(i = 0; i < n; ++ i) {
    s1_leftover[i] = lsp; s2_leftover[i] = sgroup[i].length() - (band - bh) + 1;
  } 
  // perform alignment
  for(i = 1; i <= lsp; ++ i) {   // for each iteration or query sequence position
    for(k = 0; k < n; ++ k) { // for each sequence in the target
      if(finished[k]) continue;
      rl = sgroup[k].length();
      for(x = 0; x < band; ++ x) { // for each band region in the specific target sequence
        j = i - bh + x; // computing the original index of the position
        if(j < 0) continue;
        if(j > rl) {
          if(x == 0)  {finished[k] = true;} // no need to look at the sequence any more 
          else  {
            filled[k] = true; 
            r[k] = Max3(mm[k][x - 1], mi[k][x - 1], md[k][x - 1]);
          } // record score anyway (global alignment)
          break;  // no need to continue computing for this band any more
        }
        if(j >= 1) s = f.GetMatchScore(sp[i - 1], sgroup[k][j - 1]);
        // computing order: Insertion matrix (mi) -> Match matrix (mm) -> Deletion matrix (md)
        if(j == 0)  {
          mm[k][x] = md[k][x] = md[k][x + 1] + e; mi[k][x] = mm[k][x] + g;
        } else if(x == band - 1) {
          a = mi[k][x - 1] + e; b = mm[k][x - 1] + g + e;
          mi[k][x] = a > b ? a : b;
          c = mm[k][x] + s; 
          mm[k][x] = a > c ? a : c; 
          md[k][x] = -MAX;
        } else  {
          if(x <= 0) {a = -MAX;} else {a = mi[k][x - 1] + e;}
          if(x <= 0) {b = -MAX;} else {b = mm[k][x - 1] + g + e;}
          mi[k][x] = a > b ? a : b;
          b = mm[k][x + 1] + g + e; c = md[k][x + 1] + e;
          md[k][x] = b > c ? b : c;
          t = mm[k][x] + s;
          mm[k][x] = mi[k][x] > md[k][x] ? mi[k][x] : md[k][x];
          mm[k][x] = t > mm[k][x] ? t : mm[k][x]; 
        }
      }
      if(!finished[k]) {
        -- s1_leftover[k]; -- s2_leftover[k];
      }
    }    
  }
  // check for the last time for scores that have not been filled
  // if not filled, fill it with the score that has been latestly computed
  for(k = 0; k < n; ++ k) {
    if(!filled[k]) r[k] = Max3(mm[k][band - 1], mi[k][band - 1], md[k][band - 1]);
    if(s1_leftover[k] > 0) r[k] += g + s1_leftover[k] * e;
    if(s2_leftover[k] > 0) r[k] += g + s2_leftover[k] * e;
  }
  // collect memory
  for(i = 0; i < n; ++ i) {delete [] mi[i]; delete [] mm[i]; delete [] md[i];}
  delete [] mi; delete [] mm; delete [] md;
  delete [] filled; delete [] finished;
  delete [] s1_leftover; delete [] s2_leftover;
  return;
}

int AlignBatch::AlignGlobalPairwise(
    const std::string &seq1, const std::string &seq2, 
    const int band, ScoringProt &f
) {
  int g = f.gap_open_, e = f.gap_ext_;
  if(seq1.length() == 0 && seq2.length() == 0)  return 0;
  else if(seq1.length() == 0) return g + seq2.length() * e;
  else if(seq2.length() == 0)  return g + seq1.length() * e;
  int e_band = band;
  if(e_band < 2)  e_band = 2;
  int i, j, x, a, b, c, s, t, rl = seq2.length();  
  // mm: the dp matching matrix
  int *mm = new int [band], *mi = new int [band], *md = new int [band];
  // initalizing the scores
  int bh = band / 2;
  for(i = 0; i < bh; ++ i) mm[i] = mi[i] = md[i] = -MAX;
  mm[bh] = 0; mi[bh] = md[bh] = g;
  for(i = bh + 1; i < band; ++ i)  {
    mm[i] = mi[i] = mi[i - 1] + e;
    md[i] = mm[i] + g;
  }
  bool finished = false, filled = false; int score = 0;
  
  // perform alignment
  int s1_leftover = seq1.length(), s2_leftover = seq2.length() - (band - bh) + 1;
  bool is_end_ins = false, is_end_del = false;
  for(i = 1; i <= seq1.length(); ++ i) {   // for each iteration or query sequence position
    if(finished) {
      //cout << "Jumped" << endl; 
      continue;
    } 
    
    for(x = 0; x < band; ++ x) { // for each band region in the specific target sequence
      j = i - bh + x; // computing the original index of the position
      if(j < 0) continue;
      if(j > rl) {
        if(x == 0)  {finished = true;} // no need to look at the sequence any more 
        else  {
          filled = true; 
          score = Max3(mm[x - 1], mi[x - 1], md[x - 1]);
          is_end_ins = mi[x - 1] >= score ? true : false;
          is_end_del = md[x - 1] >= score ? true : false;
        } // record score anyway (global alignment)
        break;  // no need to continue computing for this band any more
      }
      if(j >= 1) s = f.GetMatchScore(seq1[i - 1], seq2[j - 1]);
      // computing order: Insertion matrix (mi) -> Match matrix (mm) -> Deletion matrix (md)
      if(j == 0)  {
        mm[x] = md[x] = md[x + 1] + e; mi[x] = md[x] + g;
      } else if(x == band - 1) {
        md[x] = -MAX;
        a = mi[x - 1] + e;  b = mm[x - 1] + g + e;
        mi[x] = a > b ? a : b;
        c = mm[x] + s;
        mm[x] = mi[x] > c ? mi[x] : c; 
      } else  {        
        if(x <= 0) {a = -MAX;} else {a = mi[x - 1] + e;}
        if(x <= 0) {b = -MAX;} else {b = mm[x - 1] + g + e;}
        mi[x] = a > b ? a : b;
        b = mm[x + 1] + g + e; c = md[x + 1] + e;
        md[x] = b > c ? b : c;
        t = mm[x] + s;
        mm[x] = mi[x] > md[x] ? mi[x] : md[x];
        mm[x] = t > mm[x] ? t : mm[x]; 
      }
    }
    
    if(!finished) {
      -- s1_leftover; -- s2_leftover;
    }
    
    /*
    cout << "========== loop" << endl;
    for(int ix = 0; ix < band; ++ ix) {
      cout << mm[ix] << "\t";
    }
    cout << endl;
    for(int ix = 0; ix < band; ++ ix) {
      cout << mi[ix] << "\t";
    }
    cout << endl;
    for(int ix = 0; ix < band; ++ ix) {
      cout << md[ix] << "\t";
    }
    cout << endl;
    */
    
  }
  
  //cout << "left over: " << s1_leftover << "  " << s2_leftover << endl;
  //cout << "check end gap: " << is_end_ins << "  " << is_end_del << endl;
  // check for the last time for scores that have not been filled
  // if not filled, fill it with the score that has been latestly computed
  if(!filled) score = Max3(mm[band - 1], mi[band - 1], md[band - 1]);
  if(s1_leftover > 0) score += g + s1_leftover * e;
  if(s2_leftover > 0) score += g + s2_leftover * e;
  // collect memory
  delete [] mi; delete [] mm; delete [] md;
  return score;
}


int AlignBatch::AlignGlobalSingle(
    const std::string &seq1, const std::string &seq2, 
    ScoringProt &f, std::pair<std::string, std::string> &alignment
) {
  // check boundary conditions
  int g = f.gap_open_, e = f.gap_ext_;
  if(seq1.length() == 0 && seq2.length() == 0)  {
    return 0;
  } else if(seq1.length() == 0) {
    alignment = std::make_pair(string(seq2.length(), '-'), seq2);
    return g + seq2.length() * e;
  } else if(seq2.length() == 0) {
    alignment = std::make_pair(seq1, string(seq1.length(), '-'));
    return g + seq1.length() * e;
  }
  // initialize the matrices
  // mm: matix matching; mt: matrix traceback; mi: matrix insertion; md: matrix deletion
  int i, j, k, a, b, c, s, t;
  int sl = seq1.length(), tl = seq2.length();
  int **mm = new int* [sl + 1]; int **mi = new int* [sl + 1]; int **md = new int* [sl + 1];
  int **tm = new int* [sl + 1]; int **ti = new int* [sl + 1]; int **td = new int* [sl + 1];
  for(i = 0; i <= sl; ++ i) {
    mm[i] = new int [tl + 1]; mi[i] = new int [tl + 1]; md[i] = new int [tl + 1];
    tm[i] = new int [tl + 1]; ti[i] = new int [tl + 1]; td[i] = new int [tl + 1]; 
  }
  // fill-up the matrix and compute alignment
  mm[0][0] = 0; mi[0][0] = md[0][0] = g; tm[0][0] = ti[0][0] = td[0][0] = 0;
  for(j = 1; j <= tl; ++ j)  {
    mm[0][j] = mi[0][j] = mi[0][j - 1] + e; md[0][j] = mm[0][j] + g; 
    tm[0][j] = ti[0][j] = td[0][j] = 0;
  }
  for(i = 1; i <= sl; ++ i) {
    for(j = 0; j <= tl; ++ j) {
      if(j == 0)  {
        mm[i][j] = md[i][j] = md[i - 1][j] + e; mi[i][j] = mm[i][j] + g;
        tm[i][j] = ti[i][j] = td[i][j] = 0;
      } else  {        
        // insertion matrix       
        a = mi[i][j - 1] + e; b = mm[i][j - 1] + g + e;
        mi[i][j] = a > b ? a : b; ti[i][j] = a > b ? 1 : 2;
        // deleteion matrix
        b = md[i - 1][j] + e; c = mm[i - 1][j] + g + e;
        md[i][j] = b > c ? b : c; td[i][j] = b > c ? 3 : 2;
        // matching matrix
        mm[i][j] = Max3(
            mm[i - 1][j - 1] + f.GetMatchScore(seq1[i - 1], seq2[j - 1]), 
            mi[i][j], md[i][j]
        );
        // trace-back matrix
        if(mi[i][j] >= mm[i][j])  {  tm[i][j] = 1; }  // insertion
        else if(md[i][j] >= mm[i][j]) { tm[i][j] = 3;  }  // deletion 
        else {  tm[i][j] = 2; }  // matching      
      }
    }
  }
  int init_m = tm[sl][tl];
  int score = mm[sl][tl];
  // trace-back
  i = sl; j = tl; 
  string aln_seq1 = ""; string aln_seq2 = "";
  int indicator;
  switch(init_m)  {
    case 1: indicator = ti[i][j]; break;
    case 2: indicator = tm[i][j]; break;
    case 3: indicator = td[i][j]; break;
    default: cout << "Error:: AlignBatch::AlignGlobalSingle: Unrecognized trace-back initialization matrix, abort." << endl; exit(1);
  }
  while(indicator != 0 && i > 0 && j > 0) {
    //cout << "init_m:  " << init_m << endl;
    if(init_m == 1)  {
      indicator = ti[i][j];
      aln_seq1 += '-'; aln_seq2 += seq2[j - 1]; -- j;
    } else if(init_m == 2) {
      indicator = tm[i][j];
      switch(indicator)  {
        case 1: // insertion
          aln_seq1 += '-'; aln_seq2 += seq2[j - 1]; -- j; indicator = ti[i][j]; break;
        case 2: // matching
          aln_seq1 += seq1[i - 1]; aln_seq2 += seq2[j - 1]; -- i; -- j; indicator = tm[i][j]; break;
        case 3: // deletion
          aln_seq1 += seq1[i - 1]; aln_seq2 += '-'; -- i; indicator = td[i][j]; break;
        default:
          cout << "Error: AlignBatch::AlignGlobalSingle: Unrecognized trace-back matrix label, abort" << endl;
          exit(1);
      }
    } else if(init_m == 3) {
      indicator = td[i][j];
      aln_seq1 += seq1[i - 1]; aln_seq2 += '-'; -- i;
    } else  {
      cout << "Error: AlignBatch::AlignGlobalSingle: Unrecognized trace-back matrix label, abort" << endl;
      exit(1);
    }
    init_m = indicator;
  }
  
  //cout << "alignment cut: " << endl;
  
  if(i > 0)  {
    for(i; i > 0; -- i) { aln_seq1 += seq1[i - 1]; aln_seq2 += '-';  }
  } else if(j > 0) {
    for(j; j > 0; -- j) { aln_seq1 += '-'; aln_seq2 += seq2[j - 1];  }
  }

  for(auto it = aln_seq1.rbegin(); it != aln_seq1.rend(); ++ it) {  alignment.first += *it; }
  for(auto it = aln_seq2.rbegin(); it != aln_seq2.rend(); ++ it) {  alignment.second += *it; }  
  // collect the memory
  for(i = 0; i <= sl; ++ i) {
    delete [] mm[i]; delete [] mi[i]; delete [] md[i];
    delete [] tm[i]; delete [] ti[i]; delete [] td[i];
  }
  delete [] mm; delete [] mi; delete [] md;
  delete [] tm; delete [] ti; delete [] td; 
  return score;
}


/****************************************************************/
/*********************Local Alignment Modules********************/
/****************************************************************/

void AlignBatch::MultiAlignLocal(
    const int threads, char *sp, const int n, char **sgroup, 
    const int band, ScoringProt &f, std::vector<int> &r
) {
  int i, j, a, b, l = strlen(sp);
  // define the effective band: if the sequence in the target is too long
  // we cannot guarantee that the banding is enough for covering the entire target sequence
  // need to adjust it
  vector<int> eff_band(threads, 0);
  // setting the effective length through refering to the longest sequence int the set
  vector<int> range;
  range.push_back(-1);
  int chunk_size = n / threads;
  for(i = 0; i < threads - 1; ++ i) range.push_back(range.back() + chunk_size);
  range.push_back(n - 1);
  // finding the longest length of the chunk and recompute the effective banding
  for(i = 0; i < threads; ++ i) {
    for(j = range[i] + 1; j <= range[i + 1]; ++ j) {
      // now eff_band[i] stores the longest sequence length in the ith chuck
      a = abs(l - (int) strlen(sgroup[j]));
      eff_band[i] = eff_band[i] > a ? eff_band[i] : a;
    }
    eff_band[i] *= 2; eff_band[i] += band;
  } 
  // execute the alignment in multi-threaded mode
  vector<vector<int> > r_multi(threads);
  #pragma omp parallel num_threads(threads)
  {
    #pragma omp for
    for(i = 0; i < threads; ++ i) {
      AlignLocal(
          sp, range[i + 1] - range[i], &sgroup[range[i] + 1], 
          eff_band[i], f, r_multi[i]
      );
    }
    #pragma omp taskwait
  }
  // merge the alignment score results
  r.resize(n); a = 0;
  for(i = 0; i < threads; ++ i) {
    for(j = 0; j < r_multi[i].size(); ++ j) r[a ++] = r_multi[i][j];
  }
  return;
}

void AlignBatch::MultiAlignLocal(
    const int threads, std::string &sp, const int n, std::vector<std::string> &sgroup, 
    const int band, ScoringProt &f, std::vector<int> &r
) {
  int i, j, a, b, l = sp.length();
  // define the effective band: if the sequence in the target is too long
  // we cannot guarantee that the banding is enough for covering the entire target sequence
  // need to adjust it
  vector<int> eff_band(threads, 0);
  // setting the effective length through refering to the longest sequence int the set
  vector<int> range;
  range.push_back(-1);
  int chunk_size = n / threads;
  for(i = 0; i < threads - 1; ++ i) range.push_back(range.back() + chunk_size);
  range.push_back(n - 1);
  // finding the longest length of the chunk and recompute the effective banding
  for(i = 0; i < threads; ++ i) {
    for(j = range[i] + 1; j <= range[i + 1]; ++ j) {
      // now eff_band[i] stores the longest sequence length in the ith chuck
      a = abs(l - (int) sgroup[j].length());
      eff_band[i] = eff_band[i] > a ? eff_band[i] : a;
    }
    eff_band[i] *= 2; eff_band[i] += band;
  } 
  // execute the alignment in multi-threaded mode
  vector<vector<int> > r_multi(threads);
  #pragma omp parallel num_threads(threads)
  {
    #pragma omp for
    for(i = 0; i < threads; ++ i) {
      vector<string> schunk(sgroup.begin() + range[i] + 1, sgroup.begin() + range[i + 1] + 1);
      AlignLocal(
          sp, range[i + 1] - range[i], schunk, 
          eff_band[i], f, r_multi[i]
      );
    }
    #pragma omp taskwait
  }
  // merge the alignment score results
  r.resize(n); a = 0;
  for(i = 0; i < threads; ++ i) {
    for(j = 0; j < r_multi[i].size(); ++ j) r[a ++] = r_multi[i][j];
  }
  return;
}

void AlignBatch::MultiAlignLocalFull(
    const int threads, std::string &sp, const int n, 
    std::vector<std::string> &sgroup, ScoringProt &f, std::vector<FullAlnType> &r
) {
  int i; r.resize(sgroup.size());
  // execute the alignment in multi-threaded mode
  #pragma omp parallel num_threads(threads)
  {
    #pragma omp for
    for(i = 0; i < sgroup.size(); ++ i) {
      r[i].score = AlignLocalSingle(
          sp, sgroup[i], f, 
          r[i].alignment, r[i].q_interval, r[i].t_interval
      );
      r[i].sequences.first = sp; r[i].sequences.second = sgroup[i];
    }
    #pragma omp taskwait
  }
  // sort the results based on alignment score
  sort(r.begin(), r.end(), SortFullAlnType);
  return;
}


void AlignBatch::MultiAlignLocalPairwise(
    const int threads, 
    std::vector<std::string> &sgroup1, std::vector<std::string> &sgroup2,
    const int band, ScoringProt &f,
    std::vector<int> &r
) {
  if(sgroup1.size() != sgroup2.size())  {
    cout << "Warning: AlignBatch::MultiAlignLocalPairwise: number of sequences in each group do not match!!!, abort." << endl;
    exit(1);
  }
  int i, j, n = sgroup1.size();
  // computing the effective band size
  vector<int> eff_band(n, band);
  for(i = 0; i < n; ++ i) {
    eff_band[i] = abs((int) sgroup1[i].size() - (int) sgroup2[i].size()) * 2 + band;
  }
  r.resize(n);
  // compute the alignments
  #pragma omp parallel num_threads(threads)
  {
    #pragma omp for
    for(i = 0; i < n; ++ i) {
      r[i] = AlignLocalPairwise(sgroup1[i], sgroup2[i], band, f);
    }
    #pragma omp taskwait
  }
  return;
}


void AlignBatch::AlignLocal(
    char *sp, const int n, char **sgroup, 
    const int band, ScoringProt &f, 
    std::vector<int> &r
) {
  int lsp = strlen(sp);
  if(lsp <= 0 || n <= 0 || band < 2)  return;
  int i, j, k, x, t, a, b, c, s = -MAX, rl;
  int g = f.gap_open_, e = f.gap_ext_;
  // declaring the checking matrices
  bool *filled = new bool [n]; bool *finished = new bool [n]; r.resize(n, -MAX);
  // mm: matrix match; mi: matrix insert; md: matrix delete
  int **mm = new int* [n], **mi = new int* [n], **md = new int*[n];
  for(i = 0; i < n; ++ i) {
    filled[i] = finished[i] = false;
    mm[i] = new int [band]; mi[i] = new int [band]; md[i] = new int [band];
  }
  // initalizing the scores
  int bh = band / 2;
  for(i = 0; i < bh; ++ i) mm[0][i] = mi[0][i] = md[0][i] = -MAX;
  for(i = bh; i < band; ++ i)  {
    mm[0][i] = mi[0][i] = md[0][i] = 0;
  }
  int z = sizeof(int) * band;
  for(i = 1; i < n; ++ i) {
    memcpy(mm[i], mm[0], z); memcpy(mi[i], mi[0], z); memcpy(md[i], md[0], z); 
  }  
  // perform alignment
  for(i = 1; i <= lsp; ++ i) {   // for each iteration or query sequence position
    for(k = 0; k < n; ++ k) { // for each sequence in the target
      if(finished[k]) continue;
      rl = strlen(sgroup[k]);
      for(x = 0; x < band; ++ x) { // for each band region in the specific target sequence
        j = i - bh + x; // computing the original index of the position
        if(j < 0) continue;
        if(j > rl) {
          if(x == 0)  {finished[k] = true;} // no need to look at the sequence any more 
          else  {filled[k] = true;} // record score anyway (global alignment)
          break;  // no need to continue computing for this band any more
        }
        if(j >= 1) s = f.GetMatchScore(sp[i - 1], sgroup[k][j - 1]);
        // computing order: Insertion matrix (mi) -> Match matrix (mm) -> Deletion matrix (md)
        if(j == 0)  {
          mm[k][x] = md[k][x] = mi[k][x] = 0;
        } else if(x == band - 1) {
          a = mi[k][x - 1] + e; b = mm[k][x - 1] + g + e;
          mi[k][x] = a > b ? a : b;
          c = mm[k][x] + s;
          mm[k][x] = a > c ? a : c; 
          md[k][x] = -MAX;
        } else  {
          if(x <= 0) {a = -MAX;} else {a = mi[k][x - 1] + e;}
          if(x <= 0) {b = -MAX;} else {b = mm[k][x - 1] + g + e;}
          mi[k][x] = a > b ? a : b;
          b = mm[k][x + 1] + g + e; c = md[k][x + 1] + e;
          md[k][x] = b > c ? b : c;
          t = mm[k][x] + s;
          mm[k][x] = mi[k][x] > md[k][x] ? mi[k][x] : md[k][x];
          mm[k][x] = t > mm[k][x] ? t : mm[k][x]; 
        }
        // local alignment, reset score
        int max_score = Max3(mm[k][x], mi[k][x], md[k][x]);
        if(max_score < 0) {  
          mm[k][x] = mi[k][x] = md[k][x] = 0;
        }
        // record the highest score achieved so far for local alignment
        r[k] = max_score > r[k] ? max_score : r[k]; 
      }
    }
  }
  // collect memory
  for(i = 0; i < n; ++ i) {delete [] mi[i]; delete [] mm[i]; delete [] md[i];}
  delete [] mi; delete [] mm; delete [] md;
  delete [] filled; delete [] finished;
  return;
}

void AlignBatch::AlignLocal(
    std::string &sp, const int n, std::vector<std::string> &sgroup, 
    const int band, ScoringProt &f, std::vector<int> &r
) {
  int lsp = sp.length();
  if(lsp <= 0 || n <= 0 || band < 2)  return;
  int i, j, k, x, t, a, b, c, s = -MAX, rl;
  int g = f.gap_open_, e = f.gap_ext_;
  // declaring the checking matrices
  bool *filled = new bool [n]; bool *finished = new bool [n]; r.resize(n, -MAX);
  // mm: matrix match; mi: matrix insert; md: matrix delete
  int **mm = new int* [n], **mi = new int* [n], **md = new int*[n];
  for(i = 0; i < n; ++ i) {
    filled[i] = finished[i] = false;
    mm[i] = new int [band]; mi[i] = new int [band]; md[i] = new int [band];
  }
  // initalizing the scores
  int bh = band / 2;
  for(i = 0; i < bh; ++ i) mm[0][i] = mi[0][i] = md[0][i] = -MAX;
  for(i = bh; i < band; ++ i)  {
    mm[0][i] = mi[0][i] = md[0][i] = 0;
  }
  int z = sizeof(int) * band;
  for(i = 1; i < n; ++ i) {
    memcpy(mm[i], mm[0], z); memcpy(mi[i], mi[0], z); memcpy(md[i], md[0], z); 
  }  
  // perform alignment
  for(i = 1; i <= lsp; ++ i) {   // for each iteration or query sequence position
    for(k = 0; k < n; ++ k) { // for each sequence in the target
      if(finished[k]) continue;
      rl = sgroup[k].length();
      for(x = 0; x < band; ++ x) { // for each band region in the specific target sequence
        j = i - bh + x; // computing the original index of the position
        if(j < 0) continue;
        if(j > rl) {
          if(x == 0)  {finished[k] = true;} // no need to look at the sequence any more 
          else  {filled[k] = true;} // record score anyway (global alignment)
          break;  // no need to continue computing for this band any more
        }
        if(j >= 1) s = f.GetMatchScore(sp[i - 1], sgroup[k][j - 1]);
        // computing order: Insertion matrix (mi) -> Match matrix (mm) -> Deletion matrix (md)
        if(j == 0)  {
          //mi[k][x] = mi[k][x + 1] + e; mm[k][x] = md[k][x] = md[k][x + 1] + e;
          mm[k][x] = md[k][x] = mi[k][x] = 0;
        } else if(x == band - 1) {
          a = mi[k][x - 1] + e; b = mm[k][x - 1] + g + e;
          mi[k][x] = a > b ? a : b;
          c = mm[k][x] + s;
          mm[k][x] = a > c ? a : c; 
          md[k][x] = -MAX;
        } else  {          
          if(x <= 0) {a = -MAX;} else {a = mi[k][x - 1] + e;}
          if(x <= 0) {b = -MAX;} else {b = mm[k][x - 1] + g + e;}
          mi[k][x] = a > b ? a : b;
          b = mm[k][x + 1] + g + e; c = md[k][x + 1] + e;
          md[k][x] = b > c ? b : c;
          t = mm[k][x] + s;
          mm[k][x] = mi[k][x] > md[k][x] ? mi[k][x] : md[k][x];
          mm[k][x] = t > mm[k][x] ? t : mm[k][x];        
        }
        // local alignment, reset score
        int max_score = Max3(mm[k][x], mi[k][x], md[k][x]);
        if(max_score < 0) {  
          mm[k][x] = mi[k][x] = md[k][x] = 0;
        }
        // record the highest score achieved so far for local alignment
        r[k] = max_score > r[k] ? max_score : r[k];
      }
    }
  }
  // collect memory
  for(i = 0; i < n; ++ i) {delete [] mi[i]; delete [] mm[i]; delete [] md[i];}
  delete [] mi; delete [] mm; delete [] md;
  delete [] filled; delete [] finished;
  return;
}

int AlignBatch::AlignLocalPairwise(
    const std::string &seq1, const std::string &seq2, 
    const int band, ScoringProt &f
) {
  if(seq1.length() == 0 || seq2.length() == 0)  return 0;
  int e_band = band;
  if(e_band < 2)  e_band = 2;
  int i, j, x, a, b, c, s, t, ed = MAX, rl = seq2.length();
  int g = f.gap_open_, e = f.gap_ext_;
  // mm: the dp matching matrix
  int *mm = new int [band], *mi = new int [band], *md = new int [band];
  // initalizing the scores
  int bh = band / 2;
  for(i = 0; i < bh; ++ i) mm[i] = mi[i] = md[i] = -MAX;
  for(i = bh; i < band; ++ i)  {
    mm[i] = mi[i] = md[i] = 0;
  }
  bool finished = false, filled = false; int score = 0;
  // perform alignment
  for(i = 1; i <= seq1.length(); ++ i) {   // for each iteration or query sequence position
    if(finished) continue;
    for(x = 0; x < band; ++ x) { // for each band region in the specific target sequence
      j = i - bh + x; // computing the original index of the position
      if(j < 0) continue;
      if(j > rl) {
        if(x == 0)  {finished = true;} // no need to look at the sequence any more 
        else  {filled = true;} // record score anyway (global alignment)
        break;  // no need to continue computing for this band any more
      }
      if(j >= 1) s = f.GetMatchScore(seq1[i - 1], seq2[j - 1]);
      // computing order: Insertion (mi) -> Deletion (md) -> Match (mm)
      if(j == 0)  {
        mm[x] = md[x] = mi[x] = 0;
      } else if(x == band - 1) {
        a = mi[x - 1] + e; b = mm[x - 1] + g + e;
        mi[x] = a > b ? a : b;
        c = mm[x] + s;
        mm[x] = a > c ? a : c; 
        md[x] = -MAX;
      } else  {
        if(x <= 0) {a = -MAX;} else {a = mi[x - 1] + e;}
        if(x <= 0) {b = -MAX;} else {b = mm[x - 1] + g + e;}
        mi[x] = a > b ? a : b;
        b = mm[x + 1] + g + e; c = md[x + 1] + e;
        md[x] = b > c ? b : c;
        t = mm[x] + s;
        mm[x] = mi[x] > md[x] ? mi[x] : md[x];
        mm[x] = t > mm[x] ? t : mm[x];
      }      
      // local alignment, reset score
      int max_score = Max3(mm[x], mi[x], md[x]);
      if(max_score < 0) {  
        mm[x] = mi[x] = md[x] = 0;
      }
      // record the highest score achieved so far for local alignment
      score = max_score > score ? max_score : score;
    }
  }
  // collect memory
  delete [] mi; delete [] mm; delete [] md;
  return score;
}

int AlignBatch::AlignLocalSingle(
    const std::string &seq1, const std::string &seq2, 
    ScoringProt &f, std::pair<std::string, std::string> &alignment,
    std::pair<int, int> &q_interval, std::pair<int, int> &t_interval
) {
    // check boundary conditions
  if(seq1.length() == 0 || seq2.length() == 0)  {
    q_interval = t_interval = make_pair(0, 0);
    return 0;
  }
  // initialize the matrices
  // mm: matix matching; mt: matrix traceback; mi: matrix insertion; md: matrix deletion
  int i, j, k, a, b, c, s, t;
  int g = f.gap_open_, e = f.gap_ext_;
  int sl = seq1.length(), tl = seq2.length();
  //cout << "lengths: " << sl << "  " << tl << endl;
  int **mm = new int* [sl + 1]; int **mi = new int* [sl + 1]; int **md = new int* [sl + 1];
  int **tm = new int* [sl + 1]; int **ti = new int* [sl + 1]; int **td = new int* [sl + 1];
  for(i = 0; i <= sl; ++ i) {
    mm[i] = new int [tl + 1]; mi[i] = new int [tl + 1]; md[i] = new int [tl + 1];
    tm[i] = new int [tl + 1]; ti[i] = new int [tl + 1]; td[i] = new int [tl + 1];
  }
  // fill-up the matrix and compute alignment
  int max_score = 0; 
  int max_i = 0, max_j = 0;
  for(j = 0; j <= tl; ++ j)  {
    mm[0][j] = mi[0][j] = md[0][j] = 0; tm[0][j] = ti[0][j] = td[0][j] = 0;
  }
  for(i = 1; i <= sl; ++ i) {
    for(j = 0; j <= tl; ++ j) {
      if(j == 0)  {
        mm[i][j] = mi[i][j] = md[i][j] = 0; tm[i][j] = ti[i][j] = td[i][j] = 0;
      } else  {
        // insertion matrix       
        a = mi[i][j - 1] + e; b = mm[i][j - 1] + g + e;
        mi[i][j] = a > b ? a : b; ti[i][j] = a > b ? 1 : 2;
        // deleteion matrix
        b = md[i - 1][j] + e; c = mm[i - 1][j] + g + e;
        md[i][j] = b > c ? b : c; td[i][j] = b > c ? 3 : 2;
        // matching matrix
        mm[i][j] = Max3(
            mm[i - 1][j - 1] + f.GetMatchScore(seq1[i - 1], seq2[j - 1]), 
            mi[i][j], md[i][j]
        );
        // trace-back matrix
        if(mi[i][j] >= mm[i][j])  {  tm[i][j] = 1; }  // insertion
        else if(md[i][j] >= mm[i][j]) { tm[i][j] = 3;  }  // deletion 
        else {  tm[i][j] = 2; }  // matching
      }
      if(mm[i][j] > max_score) {
        max_score = mm[i][j]; max_i = i; max_j = j;
      }  
      if(mm[i][j] < 0) {
        mm[i][j] = mi[i][j] = md[i][j] = tm[i][j] = ti[i][j] = td[i][j] = 0;
      }   
    }
  }
  int init_m = tm[max_i][max_j];
  // trace-back
  i = max_i; j = max_j; 
  string aln_seq1 = ""; string aln_seq2 = "";
  int indicator;  
  switch(init_m)  {
    case 1: indicator = ti[i][j]; break;
    case 2: indicator = tm[i][j]; break;
    case 3: indicator = td[i][j]; break;
    default: cout << "Error:: AlignBatch::AlignGlobalSingle: Unrecognized trace-back initialization matrix, abort." << endl; exit(1);
  } 
  while(indicator != 0 && i > 0 && j > 0) {
    if(init_m == 1)  {
      indicator = ti[i][j];
      aln_seq1 += '-'; aln_seq2 += seq2[j - 1]; -- j;
    } else if(init_m == 2) {
      indicator = tm[i][j];
      switch(indicator)  {
        case 1: // insertion
          aln_seq1 += '-'; aln_seq2 += seq2[j - 1]; -- j; indicator = ti[i][j]; break;
        case 2: // matching
          aln_seq1 += seq1[i - 1]; aln_seq2 += seq2[j - 1]; -- i; -- j; indicator = tm[i][j]; break;
        case 3: // deletion
          aln_seq1 += seq1[i - 1]; aln_seq2 += '-'; -- i; indicator = td[i][j]; break;
        default:
          cout << "Error: AlignBatch::AlignGlobalSingle: Unrecognized trace-back matrix label, abort" << endl;
          exit(1);
      }
    } else if(init_m == 3) {
      indicator = td[i][j];
      aln_seq1 += seq1[i - 1]; aln_seq2 += '-'; -- i;
    } else  {
      cout << "Error: AlignBatch::AlignLocalSingle: Unrecognized trace-back matrix label, abort" << endl;
      exit(1);
    }
    init_m = indicator;
  }

  for(auto it = aln_seq1.rbegin(); it != aln_seq1.rend(); ++ it) {  alignment.first += *it; }
  for(auto it = aln_seq2.rbegin(); it != aln_seq2.rend(); ++ it) {  alignment.second += *it; } 
  q_interval = make_pair(i, max_i - 1); t_interval = make_pair(j, max_j - 1);
  
  //***************************************
  /*
  cout << "matrix MM: " << endl;
  for(i = 0; i <= sl; ++ i) {
    for(j = 0; j <= tl; ++ j) {
      cout << mm[i][j] << "\t";
    }
    cout << endl;
  }
  
  
  cout << "matrix MI: " << endl;
  for(i = 0; i <= sl; ++ i) {
    for(j = 0; j <= tl; ++ j) {
      cout << mi[i][j] << "\t";
    }
    cout << endl;
  }
  
  
  cout << "matrix MD: " << endl;
  for(i = 0; i <= sl; ++ i) {
    for(j = 0; j <= tl; ++ j) {
      cout << md[i][j] << "\t";
    }
    cout << endl;
  }
  
  
  cout << "matrix TM: " << endl;
  for(i = 0; i <= sl; ++ i) {
    for(j = 0; j <= tl; ++ j) {
      cout << tm[i][j] << "\t";
    }
    cout << endl;
  }
  
  
  cout << "matrix TI: " << endl;
  for(i = 0; i <= sl; ++ i) {
    for(j = 0; j <= tl; ++ j) {
      cout << ti[i][j] << "\t";
    }
    cout << endl;
  }
  
  
  cout << "matrix TD: " << endl;
  for(i = 0; i <= sl; ++ i) {
    for(j = 0; j <= tl; ++ j) {
      cout << td[i][j] << "\t";
    }
    cout << endl;
  }
  */
  
  
  //***************************************
   
  // collect the memory
  for(i = 0; i <= sl; ++ i) {
    delete [] mm[i]; delete [] mi[i]; delete [] md[i];
    delete [] tm[i]; delete [] ti[i]; delete [] td[i];
  }
  delete [] mm; delete [] mi; delete [] md;
  delete [] tm; delete [] ti; delete [] td;
  return max_score;  
}

int AlignBatch::AlignGlobalSingleEditDist(const std::string &s, const std::string &t, int band) {
  if(s.length() == 0) return t.length();
  else if(t.length() == 0)  return s.length();
  if(band < 2)  band = 2;
  int i, j, x, a, b, c, ed = MAX, rl = t.length();
  // mm: the dp matching matrix
  int *mm = new int [band];
  // initalizing the scores
  int bh = band / 2;
  for(i = 0; i < bh; ++ i) mm[i] = MAX;
  mm[bh] = 0;
  for(i = bh + 1; i < band; ++ i)  mm[i] = mm[i - 1] + 1;
  bool finished = false, filled = false; int edit_dist = 0;
  // perform alignment
  for(i = 1; i <= s.length(); ++ i) {   // for each iteration or query sequence position
    if(finished) break;
    for(x = 0; x < band; ++ x) { // for each band region in the specific target sequence
      j = i - bh + x; // computing the original index of the position
      if(j < 0) continue;
      if(j > rl) {
        if(x == 0)  {finished = true;} // no need to look at the sequence any more 
        else  {filled = true; edit_dist = mm[x - 1];} // record score anyway (global alignment)
        break;  // no need to continue computing for this band any more
      }
      if(j >= 1) ed = (s[i - 1] == t[j - 1]) ? 0 : 1;
      // computing order: Insertion matrix (mi) -> Match matrix (mm) -> Deletion matrix (md)
      if(j == 0)  {
        mm[x] = mm[x + 1] + 1;
      } else if(j == rl - 1) {
        if(x <= 0) {a = MAX;} else {a = mm[x - 1] + 1;}
        c = mm[x] + ed; mm[x] = a < c ? a : c; 
      } else  {
        if(x <= 0) {a = MAX;} else {a = mm[x - 1] + 1;}
        b = mm[x] + ed; c = mm[x + 1] + 1;
        mm[x] = a < b ? a : b; mm[x] = mm[x] < c ? mm[x] : c;
      }
    }
  }
  // check for the last time for scores that have not been filled
  // if not filled, fill it with the score that has been latestly computed
  if(!filled) edit_dist = mm[band - 1];
  // collect memory 
  delete [] mm; 
  return edit_dist;
}

int AlignBatch::KmerSimilarity(const int n, std::string &seq1, std::string &seq2) {
  if(n > seq1.length() || n > seq2.length())  return 999999;
  int i;
  unordered_map<string, int> s1_mers;
  unordered_map<string, int> s2_mers;
  for(i = 0; i <= seq1.length() - n; ++ i) {
    string mer = seq1.substr(i, n);
    if(s1_mers.find(mer) != s1_mers.end())  s1_mers[mer] ++;
    else  s1_mers[mer] = 1;
  }
  for(i = 0; i <= seq2.length() - n; ++ i) {
    string mer = seq2.substr(i, n);
    if(s2_mers.find(mer) != s2_mers.end())  s2_mers[mer] ++;
    else  s2_mers[mer] = 1;
  }
  int dist = 0;
  for(auto it = s1_mers.begin(); it != s1_mers.end(); ++ it) {
    string mer = it->first;
    if(s2_mers.find(mer) == s2_mers.end())  {
      dist += it->second * it->second;
    } else  {
      dist += (it->second - s2_mers[mer]) * (it->second - s2_mers[mer]); 
      s2_mers.erase(mer);
    }
  }
  for(auto it = s2_mers.begin(); it != s2_mers.end(); ++ it) {
    dist += it->second * it->second;
  }
  return dist;
}

int AlignBatch::WeightedKmerSimilarity(
    const int n, ScoringProt &f, std::string &seq1, std::string &seq2
) {
  if(n > seq1.length() || n > seq2.length())  return 999999;
  int i;
  unordered_map<string, int> s1_mers;
  unordered_map<string, int> s2_mers;
  for(i = 0; i <= seq1.length() - n; ++ i) {
    string mer = seq1.substr(i, n);
    if(s1_mers.find(mer) != s1_mers.end())  s1_mers[mer] ++;
    else  s1_mers[mer] = 1;
  }
  for(i = 0; i <= seq2.length() - n; ++ i) {
    string mer = seq2.substr(i, n);
    if(s2_mers.find(mer) != s2_mers.end())  s2_mers[mer] ++;
    else  s2_mers[mer] = 1;
  }
  int dist = 0;
  for(auto it = s1_mers.begin(); it != s1_mers.end(); ++ it) {
    string mer = it->first;
    int s = f.GetMaxScore(mer);
    if(s2_mers.find(mer) == s2_mers.end())  {
      dist += s * it->second * it->second;
    } else  {
      dist += s * (it->second - s2_mers[mer]) * (it->second - s2_mers[mer]); 
      s2_mers.erase(mer);
    }
  }
  for(auto it = s2_mers.begin(); it != s2_mers.end(); ++ it) {
    int s = f.GetMaxScore(it->first);
    dist += s * it->second * it->second;
  }
  return dist;
}
