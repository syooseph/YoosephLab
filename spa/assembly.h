/*
 * assembly.h
 *
 *  Created on: Feb 15, 2012
 *      Author: Youngik Yang
 */

#ifndef ASSEMBLY_H_
#define ASSEMBLY_H_

#include <string>
#include "kmer.h"
#include "graph.h"
#include "msa.h"
#include "coverage.h"
#include "container.h"
#include "eval.h"
#include "sequence.h"
#include "filter.h"
#include "semiglobal.h"
//#include "ungraph.h"
#include "PathIndex.h"
#include "Map.h"
//#include "KmerPathTable.h"

const int NOT_OVERLAP = -1000000;
const double READ_FILTER = 0.9;
const size_t PATH_OVERLAP_REGION = 30;

/**
 * Sequence strand.
 */
enum strands { POSITIVE, /**< positive strand. */
			   NEGATIVE, /**< negative strand. */
			   NONSENSE  /**< otherwise        */
};

/**
 * Path search stop criteria.
 */
enum CONDITION { NONCOMMON,   /**< No common reads exits with d distance node */
				 PATHEND,     /**< Path ends by either sink node or stop condon */
				 EXTENDFAIL,  /**< Extend fails because all nodes are weakly supported */
				 SELFLOOP,    /**< Loops itsef */
				 WEAKCURRENT, /**< Weak depth of current node */
				 WEAKMAX,     /**< Weak depth of maximum neighbor node */
				 CYCLICPATH   /**< Cycle */
};

enum PATHJOIN { PATHMERGE,
				PATHLATCH,
				MERGELATCH 
}; 
				

/**
 * String representations of stop criteria
 */
const std::string conditions[] = { "nc", "pe", "ef", "sl", "wc", "wm", "cp"};

/**
 * Relationship between pair of path.
 */
enum LINKTYPE { SUBSTRING,  /**< Substring of the other */
				BUBBLE,     /**< Bubble: slight mismatch in middle */
				SPUR,       /**< Sequence mismatch in one end  */
				SPUR_LEFT,       /**< Sequence mismatch in one end  */
				SPUR_RIGHT,       /**< Sequence mismatch in one end  */
				WEAK_SPUR,  /**< Spur with longer ends */
				LATCH,      /**< Sequence overlap in ends of pair */
				WEAK_LATCH, /**< Latch with longer ends */
				ROPE,       /**< Mismatch in both ends     */
				WEAK_ROPE,  /**< Rope with longer ends */
				LATCH_LEFT,
				LATCH_LEFT_EASY,
				LATCH_LEFT_DIFF,
				LATCH_RIGHT_EASY,
				LATCH_RIGHT,
				LATCH_RIGHT_DIFF,
				NOTYPE,     /**< Not similar each other */
				MERGE
};

/**
 * Score summary
 */
struct ScoreSummary
{
	double sum;
	double min;
	double max;
	double avg;
	double med;
	ScoreSummary( double s, double m, double M, double a, double e )
	{
		sum = s; min = m; max = M; avg = a; med =e;
	}
};

/**
 * Read information that bridges two paths
 */
struct LatchRead
{
	ReadId rid; /**< Read ID */
	PathId pid; /**< Path ID */
	int spos;   /**< Position in a read sequence */
	int rpos;   /**< Position in a reference sequence */

	/**
	   Consuctor
	*/
	LatchRead( ReadId r, PathId p, int sp, int rp )
	{
		rid = r; pid = p; spos = sp; rpos = rp;
	}
};

struct ReadPos
{
	ReadId rid;
	int spos;
	int rpos;
	ReadPos( ReadId r, int sp, int rp )
	{
		rid = r; spos = sp; rpos = rp;
	}
};

struct PairedPath
{
	PathId pid;
	int direction;
	int overlap;
	ReadIdArray query_reads;
	ReadIdArray match_reads;
	ReadIdArray latch_reads;
	std::vector<int> inits;
	PairedPath() { pid = NOT_PATH; overlap = -1000000; }
	PairedPath( PathId p, int d, size_t o, ReadIdArray q, ReadIdArray m, ReadIdArray l, std::vector<int> i ) {
		pid = p; direction = d; overlap = o; query_reads = q, match_reads = m; latch_reads = l; inits = i;
	}
};

struct PathSize
{
	size_t size;
	PathId path;

	PathSize( size_t s, PathId p )
	{
		size = s;
		path = p;
	}

/* 	bool compare( PathSize &other ) */
/* 	{ */
/* 		return size < other.size ? true : false; */

/* 	} */
};



typedef std::set<ReadId> ReadSet;
typedef std::tr1::unordered_map<PathId, SpaPath*> PathToAlnMap;
typedef std::tr1::unordered_map<KmerId, std::set<PathId> > KmerToPathMap;
typedef std::multimap<size_t, PathId> PathLengthMap;
typedef std::tr1::unordered_map<ReadId, bool> ReadFlagMap;
typedef std::set<PathId> PathIdSet;
typedef std::vector<PathId> PathIdArray;
typedef std::pair<PathId, PathId> PathIdPair;

typedef int Size;
typedef std::pair<Size, PathId> SizePathPair;
typedef std::list<SizePathPair> SizePathPairList;
typedef std::pair<std::vector<int>, PathId> PosPathPair;
typedef std::list<PosPathPair> PosPathPairList;
typedef std::list<AlignSummary> AlnSummaryList;


typedef std::list<ReadPos> ReadPosList;
typedef std::list<PathId> PathIdList;
typedef std::set<PathId> PathIdSet;
typedef std::tr1::unordered_map<PathId, PathIdSet> PathPairMap;
typedef std::map<PathId, int> PathCountMap;
typedef std::list<int> PosList;
typedef std::tr1::unordered_map<PathId, ReadPosList> PathReadMap;

typedef std::tr1::unordered_map<PathId, ReadIdList> ReadSupportMap;

typedef std::tr1::unordered_map<PathId, KmerArray> PathKmersMap;

typedef std::pair<ReadId, int> ReadPathPos;
typedef std::list<ReadPathPos> ReadPathPosList;
typedef std::tr1::unordered_map<PathId, ReadPathPosList> ReadPathPosMap;

typedef std::list<int> IntList; 
struct Scaffold
{
	PathIdList paths;
	IntList    dists;
};
typedef std::list<Scaffold> ScaffoldList;


namespace assem {
	void proceed( DeBruijnGraph &graph,
				  InvertedIndex &iindex, 
				  PathToAlnMap &path2aln_map,
				  VertexToKmerMap &vertex_map,
				  CoverageMap &kmer_coverage,
				  BitString *bstrs,
				  char *strands,
				  ReadId *pairs,
				  int &nreads, 
				  PathId *used_reads,
				  std::set<KmerType>  &debug_kmers,
				  Param &param );

void makePathPairs( PathPairMap &pair_map,
                           PathIdList &left_paths, 
                           PathIdList &right_paths );

 void inspectPathPairs( PathPairMap &path_pairs,
						//std::tr1::unordered_map<PathId, std::list<std::pair<PathId, size_t> > > &pairedreads_counts,
						PathToAlnMap &path2aln_map,
						PathId *used_reads,
						char *strands,
						ReadId *pairs, 
						Param &param
						);
 
	void pairUpReads(ReadId*, int, Param &param);


	void trimReads(InvertedIndex &iiv, 
				   NodeSet &, 
				   VertexToKmerMap &,
				   Param &param );

	void checkReads( InvertedIndex &iiv,
					 DeBruijnGraph &graph,
					 VertexToKmerMap &vertex_map,
					 CoverageMap &kmer_covearge );

	void trimMap( CoverageMap &kmer_coverage,
				  VertexToKmerMap &vertex_map,
				  NodeSet &dropped, 
				  Param &param);

	void extractPaths( DeBruijnGraph &graph,
					   BitString *bstrs,
					   char *strands,
					   ReadId *pairs,
					   int &nreads, 
					   VertexToKmerMap &vertex_map,
					   InvertedIndex &iiv, 
					   CoverageMap &kmer_coverage,
					   PathId *used_reads,
					   std::set<KmerType> &debug_kmers,
					   PathToAlnMap &,
					   Param &param );

	void extractGreedyBestPath( PathType &best_path,
								ReadIdArray &good_rids,
								Vertex curr,
								DeBruijnGraph &graph,
								InvertedIndex &iindex,
								VertexToKmerMap &vertex_map,
								AdjacencyMap &ancestors,
								int &lstop,
								int &rstop,
								Param &param );

	void extendGreedyBestNode( PathType &best_path,
							   PathType &cyclic_path,
							   Vertex &cycle_start,
							   Vertex &cycle_trigger,
							   std::queue<NodeList> &rid_queue,
							   ReadIdArray &good_rids,
							   ReadIdArray &poor_rids,
							   ReadIdArray &ncyc_rids,
							   DeBruijnGraph &graph,
							   InvertedIndex &iindex,
							   VertexToKmerMap &vertex_map,
							   AdjacencyMap &ancestors,
							   int direction,
							   int &reason,
							   Param &param);

	void findGreedyBestNeighbor( NodeSet &vset,
								 DeBruijnGraph &graph,
								 AdjacencyMap &ancestors,
								 Vertex curr_node,
								 //Vertex back_node,
								 ReadIdArray &back_comm,
								 ReadIdArray &curr_rids,
								 ReadIdArray &good_rids,
								 Vertex &max_node,
								 ReadIdArray &max_comm,
								 size_t &sum_reads,
								 std::tr1::unordered_map<Vertex, ReadIdArray> &back_rids,
								 int direction,
								 VertexToKmerMap &vertex_map,
								 InvertedIndex &iindex,
								 Param &param
								 );

	//PathLengthMap getPathLengths( PathToAlnMap &, PathIdSet & );
	PathLengthMap getPathLengths( PathToAlnMap & );

	/**
	 * Merge paths (spur, bubble, rope) and extend long overlapping paths
	 */
	void mergePairedPaths(PathToAlnMap &path2aln_map, 
						  InvertedIndex &iindex,
						  PathId *used_reads,
						  BitString *seqs,
						  char *strands,
						  ReadId *,
						  int,
						  Param &param );
	

	/**
	 * Merge paths (spur, bubble, rope) and extend long overlapping paths
	 */
	void combinePaths(PathToAlnMap &path2aln_map, 
					  InvertedIndex &iindex,
					  PathId *used_reads,
					  BitString *seqs,
					  char *strands,
					  ReadId *,
					  int,
					  Param &param );

	void makePathIdMap( KmerToPathMap &pathid_map,
						PathToAlnMap &path2aln_map );

	void makePathKmersMap( PathKmersMap &pkmers_map,
						   PathToAlnMap &path2aln_map,
						   int node_size,
						   int word_size );

	bool comparePathSize( const PathSize &one, const PathSize &two );

	
	void findLatchablePath( size_t min_filter,
							PathId sbjct_pid,
							KmerArray sbjct_kmers,
							KmerPathIndex &path_index,
							PathIdSet &skip_paths,
							PathToAlnMap &path2aln_map, 
							PathAlignPairList &palign_list,
							ReadId *pairs,
							char *strands, 
							int nreads, 
							int direction,
							Param &param,
							size_t filter_kmer);

	void recruitPaths( PathIdSet &sucess_paths,
							PathToAlnMap &path2aln_map, 
							InvertedIndex &iindex,
							PathId *used_reads,
							BitString *bstrs,
							char *strands,
							ReadId *pairs,
							int nreads,
							Param &param );

	
	ReadIdArray getBridgeReads( PathToAlnMap &path2aln_map, 
								PathId query_path, 
								PathId sbjct_path, 
								InvertedIndex &iindex, 
								PathId *used_reads,
								int direction,
								int length, 
								Param &param );


	void makePathReadMap( ReadIdArray &comm_reads,
						  PathId pid,
						  BitString *bstrs,
						  PathReadMap &latch,
						  PathToAlnMap &path2aln_map,
						  Param &param,
						  int direction );
	
	bool latchPathsWithReads( PathId &sbjct_path, 
							  PairedPath &pair_path,
							  KmerToPathMap &pathid_map,
							  PathToAlnMap &path2aln_map, 
							  BitString *bstrs,
							  char *strands,
							  ReadId *pairs,
							  PathId *used_reads,
							  PathIdSet &merged_paths,
							  PathIdSet &success_paths,
							  InvertedIndex &iindex,
							  int kmer_size,
							  Param &param);
	
	void connectPairedPathsByBridgingReads( PathId *used_reads, 
											PathToAlnMap &path2aln_map,
											BitString* bstrs,
											char *strands,
											ReadId *pairs,
											InvertedIndex &iindex,
											Param &param  );

	void connectOverlappingPairedPaths( PathToAlnMap &path2aln_map, 
										BitString *bstrs,
										char *strands,
										ReadId *pairs, 
										PathId *used_reads,
										int nreads, 
										InvertedIndex &iindex,
										Param &param );
	
	void connectPairedReadsToPath( PathId *used_reads, 
								   char *strands, 
								   ReadId *pairs, 
								   BitString *bstrs, 
								   PathToAlnMap &path2aln_map, 
								   Param &param );
	
	bool attachReadsToPathEnd( ReadIdArray &reads,
							   PathId pid,
							   PathToAlnMap &path2aln_map,
							   KmerToPathMap &pathid_map,
							   BitString *bstrs,
							   PathId *used_reads, 
							   int direction, // latch direction
							   Param &param );
	
	ReadIdArray findPairedReads( PathId &pid,
								 PathId *used_reads, 
								 ReadId *pairs, 
								 char *strands, 
								 ReadId *reads,
								 size_t nreads,
								 Param &param,
								 char strand );

	ReadIdArray getLatchableReads( ReadIdArray &reads, 
								   BitString *bstrs,
								   PathId pid,
								   PathToAlnMap &path2aln_map,
								   KmerToPathMap &pathid_map,
								   int direction,
								   Param &param );

	ReadIdArray getLatchableReads( ReadIdArray &reads, 
								   BitString *bstrs,
								   PathId pid,
								   PathToAlnMap &path2aln_map,
								   KmerToPathMap &pathid_map,
								   //int direction,
								   Param &param );
	
	void connectLongOverlappingPaths(PathToAlnMap &path2aln_map, 
									 BitString *bstrs,
									 char *strands,
									 ReadId *pairs, 
									 PathId *used_reads,
									 int nreads, 
									 InvertedIndex &iindex,
									 Param &param );
	
	void connectPaths(PathToAlnMap &path2aln_map, 
					  BitString *bstrs,
					  char *strands,
					  ReadId *pairs, 
					  PathId *used_reads,
					  int nreads, 
					  InvertedIndex &iindex,
					  Param &param );

	void tunePaths(PathToAlnMap &path2aln_map, 
				   BitString *bstrs,
				   char *strands,
				   ReadId *pairs, 
				   PathId *used_reads,
				   int nreads, 
				   InvertedIndex &iindex,
				   Param &param );
	
	void scaffoldPaths(PathToAlnMap &path2aln_map, 
					   BitString *bstrs,
					   char *strands,
					   ReadId *pairs, 
					   PathId *used_reads,
					   Param &param );

	void splitReadsByStrands( ReadIdArray &neg_rids,
							  ReadIdArray &pos_rids,
							  PathId pid,
							  PathToAlnMap &path2aln_map,
							  char *strands );
	
	std::tr1::unordered_map<ReadId, int> getReadInitPosMap( SpaPath *path );

	ReadIdArray dropConflictReads( ReadIdArray &qreads, 
								   char *strands, 
								   char strand );
	
	int estimateGap( PathId sbjct_pid, 
					 PathId query_pid, 
					 ReadIdArray &qpreads, 
					 PathToAlnMap &path2aln_map,
					 BitString *bstrs,
					 char *strands,
					 ReadId *pairs,
					 int direction,
					 Param &param);
	
	void __scaffold( Scaffold &scaffold,
					 ReadIdArray &rids,
					 PathId sbjct_pid,
					 PathToAlnMap &path2aln_map,
					 BitString *bstrs,
					 ReadId *pairs,
					 char *strands, 
					 PathId *used_reads, 
					 PathIdSet &merged_paths,
					 Map<PathIdPair, bool>::Type &scanned_pairs,
					 Param &param,
					 int direction );

	/* 	PathIdSet connectPaths( PathToAlnMap &path2aln_map,  */
/* 							InvertedIndex &iindex, */
/* 							PathId *used_reads, */
/* 							BitString *bstrs, */
/* 							char *strands, */
/* 							ReadId *pairs, */
/* 							int nreads, */
/* 							Param &param ); */
	

	void mergeClusters ( PathBins &path_bins,
						PathToAlnMap &path2aln_map, 
						BitString *bstrs,
						 PathId *used_reads,
						 char *strands,
						 ReadId *pairs,
						 Param &param );
	
	void clusterPaths( PathToAlnMap &path2aln_map, 
					   InvertedIndex &iindex,
					   PathId *used_reads,
					   BitString *bstrs,
					   char *strands,
					   ReadId *pairs,
					   int nreads,
					   Param &param );

	void clusterPairedPaths( PathToAlnMap &path2aln_map, 
							 InvertedIndex &iindex,
							 PathId *used_reads,
							 BitString *bstrs,
							 char *strands,
							 ReadId *pairs,
							 int nreads,
							 Param &param );
	
	PathIdSet latchPaths(PathIdSet &path_ids,
						PathToAlnMap &path2aln_map, 
						InvertedIndex &iindex,
						PathId *used_reads,
						BitString *seqs,
						char *strands,
						ReadId *,
						int,
						Param &param,
						int join_type );


	PathIdSet mergePaths(PathIdSet &path_ids,
						PathToAlnMap &path2aln_map, 
						InvertedIndex &iindex,
						PathId *used_reads,
						BitString *seqs,
						char *strands,
						ReadId *,
						int,
						Param &param,
						 int join_type,
						 bool,
						 bool topdown);
	
	PathIdSet linkPaths(PathIdSet &path_ids,
						PathToAlnMap &path2aln_map, 
						InvertedIndex &iindex,
						PathId *used_reads,
						BitString *seqs,
						char *strands,
						ReadId *,
						int,
						Param &param,
						int join_type );
	
	std::pair<int, PathId> findMostSimilarPath( KmerId *kmers,
												size_t nkmer,
												KmerToPathMap &pathid_map);
	void addPath( KmerId *kmers,
				  size_t nkmer,
				  DeBruijnGraph &path_graph,
				  KmerToVertexMap &vertex_map );


	void addReadsToPath( PathId pid,
						 ReadPathPosList rposs,
						 PathToAlnMap &path2aln_map,
						 BitString *bstrs,
						 PathId *used_reads,
						 Param &param );

	void simplePathClean( ReadPathPosMap &added_reads,                            
						  PathIdSet &merged_paths, 
						  PathToAlnMap &path2aln_map,
						  PathId *used_reads,
						  BitString *bstrs,
						  char *strands, 
						  ReadId *pairs, 
						  Param &param);
	
	void cleanDirtyPaths( PathIdSet &dirty_paths,
						  PathIdSet &merged_paths,
						  PathToAlnMap &path2aln_map,
						  PathId *used_reads,
						  BitString *bstrs,
						  char *strands, 
						  ReadId *pairs, 
						  Param &param);

	void recruitReads( PathToAlnMap &path2aln_map, 
					   InvertedIndex &iindex,
					   BitString *bstrs,
					   char *strands,
					   ReadId *pairs,
					   PathId *used_reads,
					   size_t nreads,
					   Param &param );

	void extendShortOverlapPaths( PathToAlnMap &path2aln_map, 
								  BitString *bstrs,
								  char *strands,
								  ReadId *pairs, 
								  PathId *used_reads,
								  int nreads,
								  InvertedIndex &iindex,
								  Param &param );

	ReadIdList getUnusedReads( PathId *used_reads,
							   int nreads );

	/**
	 * Connect pairs of separated paths with support of reads that overlaps with the path pair.
	 */
	void extendPathsByBridgingReads( std::list<LatchRead> &latch_reads,
									 PathId *used_reads, 
									 PathToAlnMap &path2aln_map,
									 BitString* bstrs,
									 char *strands,
									 ReadId *pairs,
									 KmerToPathMap &pathid_map,
									 PathIdSet &merged_paths,
									 PathIdSet &successPaths,
									 InvertedIndex &,
									 Param &);
	

	void trim( DeBruijnGraph &graph, CoverageMap& kmer_coverage, VertexToKmerMap &vertex_map, InvertedIndex &iindex, Param &param );


	ScoreSummary computeEntropy( Profile &profile );
	ScoreSummary computeDepth( Profile &profile );

	void setStartMinCoverage( InvertedIndex &iindex, Param &param );

	void initPathSearchPriority( std::multimap<double, Vertex> &cov_map,
								 //std::tr1::unordered_map<Vertex, double> &cov_map,
								 DeBruijnGraph &graph,
								 InvertedIndex &iindex,
								 VertexToKmerMap &vertex_map,
								 AdjacencyMap &ancestors,
								 Param &param );
	
	bool validStartNode( Vertex &v, 
						 InvertedIndex &iindex, 
						 AdjacencyMap &ancestors, 
						 VertexToKmerMap &vertex_map, 
						 NodeSet &deleted_nodes, 
						 std::set<KmerType> &debug_kmers,
						 Param &param  );

	bool shortPath( PathType &max_path,
					VertexToKmerMap &vertex_map, 
					Param &param );
	
	void append( ReadIdList &reads, 
				 ReadId rids[],
				 size_t count ) ;


	void invalidate( ReadIdArray &good_rids,
					 ReadIdArray &poor_rids,
					 ReadIdArray &del_rids,
					 bool sorted );
	
	void updateQueue( std::queue<NodeList> &nq,
					  NodeList &nodes,                                      
					  ReadIdArray &poor_rids, 
					  InvertedIndex &iindex, 
					  VertexToKmerMap &vertex_map, 
					  size_t q_size );
	
	std::pair<Vertex, int> getTraceNode( PathType &best_path,
										 int direction,
										 Param &param );
	
	NodeSet findGoodNeighbors( DeBruijnGraph &graph,
							   Vertex curr_node, 
							   PathQueue queue,  // pass by copy
							   VertexToKmerMap &vertex_map,
							   int direction,
							   Param &param );

	ReadIdArray getNewSupports ( ReadIdArray &good_rids,
								 ReadIdArray &poor_rids, 
								 ReadIdArray &curr_rids );

	ReadIdArray trimUsedReads( ReadIdArray &path_rids,
							   PathId *used_reads );
	

	SpaPath* getPathReadPile( PathType &max_path,
							  VertexToKmerMap &vertex_map,
							  InvertedIndex &iindex,
							  BitString *bstrs,
							  ReadIdArray &good_rids, 
							  PathId *used_reads,
							  Param &param );

	PathId savePath( SpaPath *spath,
					 PathToAlnMap &path2aln_map );
	
	void printPathSummary( std::string consensus,
						   size_t kmer_depth,
						   KmerId start_kmer,
						   int nparents,
						   int nchildren,
						   int lstop,
						   int rstop,
						   ScoreSummary &coverage,
						   ScoreSummary &entropy,
						   int creads,
						   Param &param
						   );

	void adjustReadIds ( PathType &max_path, 
						 DeBruijnGraph &graph,
						 ReadIdArray &path_reads,
						 InvertedIndex &iindex, //std::tr1::unordered_map<KmerId, Read*> &kmer2read_map,
						 VertexToKmerMap &vertex_map, 
						 int k );
	
	void updateReadIds ( ReadIdArray &nrids,
						 KmerId &kid,
						 InvertedIndex &iindex ) ;
	
	
	void dropWeakNodes( NodeSet &deleted_nodes,
						PathType &max_path,
						DeBruijnGraph &graph,
						AdjacencyMap &ancestors, 
						InvertedIndex &iindex,
						VertexToKmerMap &vertex_map,
						Param &param );
	

	void dropBadPriorities(std::tr1::unordered_map<Vertex, double> &cov_map, NodeSet &nodes);

	void updatePriorities(PathType &max_path, 
						  std::tr1::unordered_map<Vertex, double> &cov_map,
						  DeBruijnGraph &graph,
						  InvertedIndex &iindex, 
						  VertexToKmerMap &vertex_map,
						  AdjacencyMap &ancestors,
						  NodeSet &deleted_nodes, 
						  Param &param );

	std::pair<Vertex, double> findMaxPriority( std::tr1::unordered_map<Vertex, double> &cov_map );
 
	void updatePathSearchPriority( Vertex v,
								   std::multimap<double, Vertex> &cov_map,
								   //std::tr1::unordered_map<Vertex, double> &cov_map,
								   DeBruijnGraph &graph,
								   InvertedIndex &iindex,
								   VertexToKmerMap &vertex_map,
								   AdjacencyMap &ancestors,
								   Param &param );
	
	void displayPathSearchSummary( int nstart, int delnodes, int npaths, int reads );

	int getMinKmerCount( int seqlen, Param &param  );

	void initPosFlags( std::vector<int> &pvec,
                   int qnkmer, 
                   int rnkmer,
                   KmerId *qkmers,
                   KmerId *rkmers,
							  std::map<int, std::set<int> > &collide );
	
	bool latchable( std::map<int, std::set<int> > &collide,
						   std::vector<int> &pvec,
						   int qnkmer, 
						   int rnkmer,
						   int &count,
						   int min_length,
						   int direction,
						   Param &param  );
	
	
	std::vector<int> updatePosFlags(std::vector<int> &pvec, int count, int rstart, int qstart );
	
	bool goodLatch(std::vector<int> &pvec, KmerId *qkmers, KmerId *rkmers, int qnkmer, int rnkmer, int direction, Param &param  );

	std::vector<int> makePosFlags( int qnkmer, 
										  int rnkmer, 
										  KmerId *qkmers,
										  KmerId *rkmers,
										  int min_length, Param &param);

	void setKmerMap( std::tr1::unordered_map<KmerId, int> & count_map, 
					 KmerArray &kmers );

	size_t countSameKmers( std::tr1::unordered_map<KmerId, int> &count_map,
						   KmerArray &query_kmers );
	
	size_t countSameKmers( std::string &sbjct,
                              std::string &query,
                              size_t small_k );

size_t countSameKmers( KmerToPathMap &pathid_map,
					   PathToAlnMap &path2aln_map,
					   PathId query_pid,
					   PathId sbjct_pid,
					   size_t small_k);

std::map<PathId, int>  getKmerCountMap( PathIdSet pids,
                                        KmerId *kmers,
                                        int nkmer,
                                        KmerToPathMap &pathid_map,
											   PathIdSet &merged_paths );

PosPathPairList getKmerPosList( std::map<PathId, int> &count_map,
                                KmerId *kmers,
                                int nkmer,
                                PathToAlnMap &path2aln_map,
                                int mink, 
								bool is_query,
								Param &param );
 
 PathIdSet searchOverlapCandidates(  //PathIdSet pids,
								PathId sbjct_pid,
								KmerId *skmers,
								//int snkmer,
								int srch_nkmer,
								PathToAlnMap &path2aln_map,
								KmerToPathMap &pathid_map,
								PathIdSet &merged_paths,
								int mink,
								Param &param,
								int direction );
 
 PosPathPairList searchSimilarPaths( PathIdSet pids,
									 KmerId *kmers,
									 int nkmer,
									 PathToAlnMap &path2aln_map,
									 KmerToPathMap &pathid_map,
									 PathIdSet &merged_paths,
									 int mink,
									 bool is_query,
									 Param &param );


void addPathIds( PathId pid,
                    KmerId *kmers,
                    size_t nkmer, 
						KmerToPathMap &pathid_map);

void dropPathIds( PathId pid,
                  KmerId *kmers,
                  size_t nkmer, 
						 KmerToPathMap &pathid_map);

bool mergePathPair( PathId query_pid, 
                    PathId sbjct_pid, 
                    AlignSummary &summary,
                    PathIdSet &merged_paths,
                    PathToAlnMap &path2aln_map,
                    //KmerToPathMap &pathid_map,
                    BitString *bstrs,
                    char *strands, 
                    ReadId *pairs,
                    PathId *used_reads,
                    InvertedIndex &iindex,
						   Param &param );

 bool inorder(std::vector<int> &pos);
 void makeTightBound(IntPair &block, std::vector<int> &pos_vec);
 std::vector<IntPair> makeBlocks( std::vector<int> &pos_vec );
 IntPair expandBound( IntPair bound, 
					  std::vector<int> &pos_vec,
					  Param &param);
 
 IntPair getMatchRange(std::vector<int>& pos_vec, Param &param );

void getMatchRanges( std::vector<IntPair> &rranges, 
                            std::vector<IntPair> &qranges, 
                            std::vector<int> &pos_vec,
                            PathId query_pid, 
                            PathId sbjct_pid, 
                            ReadId *pairs, 
                            int nreads, 
                            Param &param );
 IntPair getMatchRangeMax(std::vector<int>& pos_vec, Param &param );
 IntPair getMatchRangeLeft(std::vector<int>& pos_vec, Param &param );
 IntPair getMatchRangeRight(std::vector<int>& pos_vec, Param &param );

int countMatchKmers(std::vector<int>& pos_vec,
                     size_t s,
						   size_t e );


bool indeled( std::vector<int>& pos_vec,
					 IntPair &range );

bool ordered(std::vector<int>& pos_vec,
					IntPair &range );

 int getPathAlignType( AlignSummary &summary );
 int getMatchType( IntPair &qrange,
				   IntPair &rrange,
				   KmerId *qkmers,
				   KmerId *rkmers,
				   int     qnkmer,
				   int     rnkmer,
				   Param &param  );

int getLinkType(IntPair &qrange,
                IntPair &rrange,
                KmerId *qkmers,
                KmerId *rkmers,
                size_t qnkmer,
                size_t rnkmer,
                size_t latch_offset, 
                size_t spur_offset,
					   Param &param  );

int getPathPairType( std::vector<int> &pos_vec,
					 IntPair &qrange,
					 IntPair &rrange,
					 KmerId *qkmers,
					 KmerId *rkmers,
					 size_t qnkmer,
					 size_t rnkmer,
					 size_t latch_offset, 
					 size_t spur_offset,
					 Param &param  );
 
 
 bool isBubble(size_t qs,
			  size_t qe,
			  size_t ms,
			  size_t me,
			  size_t qnkmer,
			  size_t rnkmer,
			  KmerId *qkmers,
			  KmerId *rkmers,
			  size_t spur_offset,
			  int &type,
			  Param &param );
 
bool isSpurRight(size_t qs,
				 size_t qe,
				 size_t ms,
				 size_t me,
				 size_t qnkmer,
				 size_t rnkmer,
				 KmerId *qkmers,
				 KmerId *rkmers,
				 size_t spur_offset,
				 int &type,
				 Param &param );

bool isSpurLeft(size_t qs,
                      size_t qe,
				  size_t ms,
				  size_t me,
				  size_t qnkmer,
				  size_t rnkmer,
				  KmerId *qkmers,
				  KmerId *rkmers,
				  size_t spur_offset,
				  int &type,
				  Param &param );


bool isLatchRight(size_t qs,
                      size_t qe,
				  size_t ms,
				  size_t me,
				  size_t qnkmer,
				  size_t rnkmer,
				  KmerId *qkmers,
				  KmerId *rkmers,
				  size_t spur_offset,
				  //int &type,
				  Param &param );

bool isLatchLeft(size_t qs,
                      size_t qe,
				  size_t ms,
				  size_t me,
				  size_t qnkmer,
				  size_t rnkmer,
				  KmerId *qkmers,
				  KmerId *rkmers,
				  size_t spur_offset,
				  //int &type,
				  Param &param );

bool isFrayedRope(size_t qs,
                      size_t qe,
				  size_t ms,
				  size_t me,
				  size_t qnkmer,
				  size_t rnkmer,
				  KmerId *qkmers,
				  KmerId *rkmers,
				  size_t spur_offset,
				  int &type,
				  Param &param );

 bool latchableReadRight( size_t qs,
						  size_t qe,
						  size_t ms,
						  size_t me,
						  size_t qnkmer,
						  size_t rnkmer,
						  KmerId *qkmers,
						  KmerId *rkmers,
						  size_t spur_offset,
						  Param &param);
 
 bool latchableReadLeft( size_t qs,
						 size_t qe,
						 size_t ms,
						 size_t me,
						 size_t qnkmer,
						 size_t rnkmer,
						 KmerId *qkmers,
						 KmerId *rkmers,
						 size_t spur_offset,
						 Param &param);

	
AlignSummary compareBySubstitution( std::string &query,
                                    std::string &sbjct,
                                    int qnkmer,
                                    int rnkmer,
                                    int qs,
                                    int qe,
                                    int rs,
                                    int re,
										   Param &param  );

AlignSummary compareByAlignment( std::string query,
                                 std::string sbjct,
                                 int qnkmer,
                                 int rnkmer,
                                 int qs,
                                 int qe,
                                 int rs,
                                 int re,
										Param &param );


void adjustRanges( IntPair &qrange,
                   IntPair &rrange,
                   int qnkmer,
                   int rnkmer,
						  Param &param );

void adjustLatchRange( IntPair &qrange,
                   IntPair &rrange,
                   int qnkmer,
                   int rnkmer,
						  Param &param );


void Interchange(PathId &qid, 
                 PathId &mid, 
                 KmerId **qkmers, 
                 KmerId **rkmers, 
                 size_t &qnkmer, 
                 size_t &rnkmer, 
                 IntPair &qrange, 
                 IntPair &rrange, 
						Param &param );

 bool pickMergePath( PosPathPairList &pos_paths,
					 PathId &qid,
					 PathId &mid,
					 AlignSummary &summary,
					 PathToAlnMap &path2aln_map,
					 KmerToPathMap &pathid_map,
					 ReadId *pairs, 
					 size_t nreads,
					 Param &param  );

 void printGlobalAlignment( GlobalAlignPair &aln );
 void printLocalAlignment( LocalAlignPair &aln );
 void printPaths( PathToAlnMap &path2aln_map, Param &param  );
 void printPaths( PathToAlnMap &path2aln_map, Param &param  );
void dropMergedPaths( PathIdSet &merged_paths, 
							 PathToAlnMap &path2aln_map);


 std::string getSequence( SpaPath *spath, Param &param  );

void combinePaths( PathToAlnMap &path2aln_map, 
                   InvertedIndex &iindex,
                   //char **seqs,
                   PathId *used_reads,
                   BitString *bstrs,
                   char *strands,
                   ReadId *pairs,
                   int nreads,
						  Param &param );

 std::multimap<int, PathId> getSimilarPaths( PathId &query_pid, 
											 //PathIdSet &skip_pids, 
                                         PathToAlnMap &path2aln_map, 
                                         KmerToPathMap &pathid_map,
                                         PathIdSet &merged_paths,
                                         ReadId *pairs,
                                         int nreads,
												   Param &param );


bool searchMostSimilarPath( PathId &query_pid, 
                            PathId &sbjct_pid, 
                            AlignSummary &summary,
                            PathToAlnMap &path2aln_map, 
                            KmerToPathMap &pathid_map,
                            PathIdSet &merged_paths,
                            ReadId *pairs,
                            int nreads,
								   Param &param );


 bool alreadyCompared( PathId query_pid,
					   PathId sbjct_pid, 
					   std::tr1::unordered_map<PathId, PathIdList> &compared_paths );
 
 bool mergeSimilarPath( PathId &query_pid, 
						//PathIdSet &skip_pids,
						PathToAlnMap &path2aln_map, 
						KmerToPathMap &pathid_map,
						PathIdSet &merged_paths,
						std::tr1::unordered_map<PathId, PathIdList> &compared_paths,
						std::tr1::unordered_map<PathId, bool> &updated_paths,
						ReadId *pairs,
						char *strands, 
						int nreads,
						BitString *bstrs,
						PathId *used_reads,
						InvertedIndex &iindex,
						Param &param );


 std::pair<int, AlignSummary> getMaxAlign(PathId query_pid, 
						 PathId sbjct_pid, 
						 intVec &types, 
						 std::vector<IntPair> &rranges, 
						 std::vector<IntPair> &qranges, 
						 std::vector<int> &pos_vec,
						 PathToAlnMap &path2aln_map, 
						 Param &param);

AlignSummary alignPathPair( PathId query_pid, 
							PathId sbjct_pid, 
							int type,
							IntPair &qrange, 
							IntPair &rrange,
							std::vector<int> &pos_vec,
							PathToAlnMap &path2aln_map,
							Param &param );
 
 bool enoughKmerShare( int qnkmer, 
					   IntPair rrange,
					   std::vector<int> &pos_vec,
					   int type, 
					   Param &param );

 PosPathPairList getMergiblePaths( PathId &query_pid, 
								   PathToAlnMap &path2aln_map, 
								   KmerToPathMap &pathid_map,
								   PathIdSet &merged_paths,
								   ReadId *pairs,
								   int nreads,
								   Param &param );
 
void determineMatchRangesAndTypes( intVec &types, 
                                   std::vector<IntPair> &rranges, 
                                   std::vector<IntPair> &qranges, 
                                   std::vector<int> &pos_vec,
                                   PathId query_pid, 
                                   PathId sbjct_pid, 
                                   ReadId *pairs, 
                                   int nreads, 
                                   PathToAlnMap &path2aln_map,
										  Param &param );

ReadIdArray getCommonReads(PathId query_pid, 
						   PathId  sbjct_pid,
						   PathToAlnMap &path2aln_map, 
						   ReadId *pairs,
                           int nreads, 
						   Param &param );

intVec remakePositionVector( PathId qpid,
                             PathId spid, 
                             PathToAlnMap &path2aln_map, 
									Param &param );

bool __mergeBestPath( PathId query_pid, 
                             PathId sbjct_pid, 
                             intVec &types, 
                             std::vector<IntPair> &rranges, 
                             std::vector<IntPair> &qranges, 
                             std::vector<int> &pos_vec,
                             ReadIdArray &sbjpreads, 
                             int direction, 
                             PathToAlnMap &path2aln_map, 
                             KmerToPathMap &pathid_map,
                             BitString *bstrs,
                             char *strands, 
                             ReadId *pairs,
                             PathId *used_reads,
                             InvertedIndex &iindex,
                             PathIdSet &merged_paths,
                             Param &param);

void getMatchTypes(intVec &types, 
                   std::vector<IntPair> rranges, 
                   std::vector<IntPair> qranges, 
                   std::vector<int> &pos_vec,
                   PathId query_pid, 
                   PathId sbjct_pid, 
                   PathToAlnMap &path2aln_map,
						  Param &param );

 bool alignMiddle( AlignSummary &summary,
				   std::string sbjct,
				   std::string query,
				   //int same_kmer, 
				   int &type, 
				   Param &param );


 bool alignEnd( AlignSummary &summary,
				std::string sbjct,
				std::string query,
				PathId query_pid,
				PathId sbjct_pid,
				std::set<IntPair> &bad_range,
				PathToAlnMap &path2aln_map, 
				Param &param, 
				int anchor );
 

bool checkSimpleLatch( AlignSummary &summary, 
                       KmerId *qkmers, 
                       KmerId *rkmers, 
                       int qnkmer,
                       int rnkmer, 
                       IntPair &qrange, 
                       IntPair &rrange, 
							  Param   &param );

bool checkLatchAlign( AlignSummary &summary, 
                      std::string &query,
                      std::string &sbjct,
                      int anchor,
							 Param   &param );

 std::pair<bool, int>  __tryMergePathPair( PathId &query_pid, 
											   PathId  sbjct_pid,
										   //int same_nkmer,
											   PathToAlnMap &path2aln_map, 
										   //KmerToPathMap &pathid_map,
											   PathIdSet &merged_paths,
											   ReadId *pairs,
											   char *strands,
											   int nreads,
											   BitString *bstrs,
											   PathId *used_reads,
											   InvertedIndex &iindex,
										   Param &param,
										   int join_type);
 
/*  std::pair<bool, int> __tryMergePath( PathId &query_pid, */
/*                             PathId  sbjct_pid, */
/*                             intVec pos_vec, */
/*                             PathToAlnMap &path2aln_map, */
/*                             KmerToPathMap &pathid_map, */
/*                             PathIdSet &merged_paths, */
/*                             ReadId *pairs, */
/*                             char *strands, */
/*                             int nreads, */
/*                             BitString *bstrs, */
/*                             PathId *used_reads, */
/*                             InvertedIndex &iindex, */
/*                             Param &param ); */

 void linkPaths(PathIdSet &path_ids,
				PathToAlnMap &path2aln_map, 
				InvertedIndex &iindex,
				PathId *used_reads,
				BitString *bstrs,
				char *strands,
				ReadId *pairs,
				int nreads,
				Param &param );
 
 int getUnusedReadCount( PathId *used_reads,
						 int nreads );

 bool latchableType(int type);
 bool mergibleType(int type);
 IntPair __extractSbjectRange(std::vector<int> &pos_vec,
							  Param &param,
							  int anchor );
 bool __validRange(IntPair &rrange,
				   std::vector<int> &pos_vec, 
				   Param &param);

bool __alignReadToPath(int &type,
					   AlignSummary &summary,
                              IntPair &qrange,
                              IntPair &rrange,
                              KmerId *qkmers,
                              size_t qnkmers,
                              KmerId *skmers,
                              size_t snkmers,
                              Param &param );

 void __mergeToPath( ReadId curr_rid,
					 PathId *used_reads,
					 PathToAlnMap &path2aln_map,
					 PathId sbjct_pid,
					 IntPair &qrange,
					 IntPair &rrange,
					 AlignSummary &summary,
					 int &merged,
					 //PathIdSet &dirty_paths,
					 ReadPathPosMap &added_reads,
					 Param &param);
 
	 
void recruit( BitString *bstrs,
              PathId *used_reads,
              ReadId urid,
              PathToAlnMap &path2aln_map,
              KmerToPathMap &pathid_map,
              std::list<LatchRead> &latch_reads,
              int &merged,
			  //PathIdSet &dirty_paths,
			  ReadPathPosMap &added_reads,
			  Param &param
					 );

bool __goodReadToJoin( std::vector<int> &pos_vec,
                              IntPair &qrange,
                              IntPair &rrange,
                              int &type,
                              AlignSummary &summary,
                              KmerId *qkmers,
                              size_t qnkmers,
                              KmerId *skmers,
                              size_t snkmers,
                              Param &param,
                              int anchor );


void recruitReads( PathToAlnMap &path2aln_map,
                   InvertedIndex &iindex,
                   BitString *bstrs,
                   char *strands,
                   ReadId *pairs,
                   PathId *used_reads,
                   size_t nreads,
						  Param &param );

ReadIdList getUnusedReads( PathId *used_reads,
								  int nreads );
 ReadIdArray getReadIds( std::list<ReadPos> &rpos_list );
 std::tr1::unordered_map<ReadId, int> getReadStartPos( ReadPosList rp, int which );
void pairUpPaths( PathPairMap &pair_map,
                  PathReadMap &latch_left,
                  PathReadMap &latch_right,
                  std::list<LatchRead> &latch_reads,
						 Param &param );
void commonReadSet( PathId left,
                    PathPairMap &pair_map,
                    PathReadMap &latch_left,
                    PathReadMap &latch_right,
                    std::map<int, PathId> &commpid_map, 
                    std::tr1::unordered_map<PathId, ReadIdArray> &commset_map,
						   Param &param );

void validLatchReads( PathId *used_reads,
                      BitString *bstrs,
                      std::string &lstr,
                      std::string &rstr,
                      PathId left,
                      PathId right,
                      PathReadMap &latch_left,
                      PathReadMap &latch_right,
                      ReadIdArray &comm_reads,
                      ReadIdArray &good_reads,
                      std::vector<int> &read_inits,
                      int &overlap,
					  Param &param
							 );

void inspectLatchReads( PathId left,
                               PathId curr,
                               ReadIdArray &good_reads,
                               std::vector<int> &read_inits,
                               int &overlap, 
                               std::map<int, PathId> &commpid_map,
                               std::tr1::unordered_map<PathId, ReadIdArray> &commset_map,
                               BitString             *bstrs,
                               PathReadMap           &latch_left,
                               PathReadMap           &latch_right,
                               PathToAlnMap          &path2aln_map,
                               PathId *used_reads,                                
                               Param                 &param );

void findMaxExtendiblePath( PathId                 left,
                                   std::map<int, PathId> &commpid_map,
                               std::tr1::unordered_map<PathId, ReadIdArray> &commset_map,
                                   PathId                &mpath,
                                   ReadIdArray           &mrids,
                                   std::vector<int>      &mbegs,
                                   int                   &mover,
                                   PathId                *used_reads,
                                   BitString             *bstrs,
                                   PathReadMap           &latch_left,
                                   PathReadMap           &latch_right,
                                   PathToAlnMap          &path2aln_map,
                                   Param                 &param
                                   );

 bool connectPathPair( BitString *bstrs,
					   char *strands,
					   ReadId *pairs,
					   PathToAlnMap &path2aln_map,
					   KmerToPathMap &pathid_map,
					   PathIdSet &merged_paths,
					   PathIdSet &success_paths,
					   PathId sbjct_pid,
					   PathId query_pid,
					   PathId *used_reads,
					   ReadIdArray &good_reads,
					   std::vector<int> &read_inits,
					   int overlap,
					   int direction,
					   InvertedIndex &iindex,
					   Param &param
					   );
 
bool latch( BitString *bstrs,
            char *strands,
            ReadId *pairs,
            PathToAlnMap &path2aln_map,
            KmerToPathMap &pathid_map,
            PathIdSet &merged_paths,
            PathIdSet &success_paths,
            PathId left,
            PathId right,
            PathId *used_reads,
            ReadIdArray &good_reads,
            std::vector<int> &read_inits,
            int overlap,
            InvertedIndex &iindex,
			Param &param
				   );

void remapLatchReads( PathId left, 
                  PathId right, 
                  PathPairMap &pair_map,
                  PathReadMap &latch_left,
                  PathReadMap &latch_right,
                  ReadIdArray &good_reads,
                  int lstr_size,
						 int &overlap );

 bool joinReadsLeftOver( PathId *used_reads,
						 BitString *bstrs,
						 ReadIdArray inspected,
						 ReadIdArray path_reads, /* sorted already */
						 PathId pid,
						 PathToAlnMap &path2aln_map,
						 Param &param);
 
 void getPairedReadCounts( std::multimap<int, PathId> &pcount_map,
                           PathId left,
                           PathIdSet &pids,
                           ReadId *pairs,
                           PathId *used_reads,
						   char *strands,
						   PathToAlnMap &path2aln_map,
						   PathIdSet &merged_paths, 
						   Param &param );

void latchPair( BitString *bstrs,
                char *strands,
                ReadId *pairs,
                PathId *used_reads,
                PathPairMap &pair_map,
				//std::tr1::unordered_map<PathId, std::list<std::pair<PathId, size_t > > > &pairedreads_counts,
                PathReadMap &latch_left,
                PathReadMap &latch_right,
                PathToAlnMap &path2aln_map,
                KmerToPathMap &pathid_map,
                PathIdSet &merged_paths,
                PathId left,
                InvertedIndex &iindex,
                PathIdSet &successPaths,
				PathIdSet &dirty_paths,
					   Param &param );

void latchPairs(std::list<LatchRead> &latch_reads,
                PathId *used_reads, 
                PathToAlnMap &path2aln_map,
                BitString* bstrs,
                char *strands,
                ReadId *pairs,
                KmerToPathMap &pathid_map,
                PathIdSet &merged_paths,
                PathIdSet &successPaths,
				PathIdSet &dirty_paths,
                InvertedIndex &iindex,
                PathReadMap &latch_left, 
                PathReadMap &latch_right,
                PathPairMap &pair_map,
				//std::tr1::unordered_map<PathId, std::list<std::pair<PathId, size_t> > > &pairedreads_counts,
					   Param &param );
bool recruitToPath( ReadId rid,
                    PathId path,
                    PathId *used_reads, 
                    BitString* bstrs,
                    PathToAlnMap &path2aln_map,
                    KmerToPathMap &pathid_map,
                    PathIdSet &merged_paths,
                    PathIdSet &successPaths,
					Param &param
						   );

int recruitReadsToPath( ReadIdSet &rset,
                        PathId path,
                        PathId *used_reads, 
                        BitString* bstrs,
                        PathToAlnMap &path2aln_map,
                        KmerToPathMap &pathid_map,
                        PathIdSet &merged_paths,
                        PathIdSet &successPaths,
						Param &param
							   );

void latchReadsToPath( PathReadMap &latch_reads,
                       PathId *used_reads, 
                       char *strands,
                       ReadId *pairs, 
                       BitString* bstrs,
                       PathToAlnMap &path2aln_map,
                       KmerToPathMap &pathid_map,
                       PathIdSet &merged_paths,
                       PathIdSet &successPaths,
                       int direction,
							  Param &param );

void extendPathsByBridgingReads( std::list<LatchRead> &latch_reads,
                        PathId *used_reads, 
                        PathToAlnMap &path2aln_map,
                        BitString* bstrs,
                        char *strands,
                        ReadId *pairs,
                        KmerToPathMap &pathid_map,
                        PathIdSet &merged_paths,
                        PathIdSet &successPaths,
								 PathIdSet &dirty_paths,
								 InvertedIndex &iindex,
										Param &param  );

PathIdArray findPairedPaths( PathId query_path,
                            PathId *used_reads,
                            ReadId *pairs, 
                            ReadId *reads, 
							 size_t nreads,
									Param &param );


 ReadIdArray getPairedReads( ReadIdArray &this_reads,
							PathId &this_path,
							 PathId &that_path,
							 PathId *used_reads,
							 ReadId *pairs, 
							 PathToAlnMap &path2aln_map );
 
void setPairedReads ( ReadIdArray &qreads,
                      ReadIdArray &sreads,
                      PathId query_path,
                      PathId sbjct_path,
                      PathId *used_reads,
                      ReadId *pairs,
							 PathToAlnMap &path2aln_map );
IntPair determineStrands( ReadIdArray &preads,
                          ReadId *pairs,
								 char *strands, Param &param );

int getDirection(ReadIdArray &preads,
                 ReadId *pairs,
						char *strands, Param &param );

 IntPair getSimpleLatchRightRange( std::vector<int>& pvec );
 IntPair getSimpleLatchLeftRange( std::vector<int>& pvec, int qnkmer );
 
 IntPair checkLongOverlap(PathId qpid, PathId spid, PathToAlnMap &path2aln_map, int min_length, int direction, Param &param  );

IntPair checkShortOverlap( PathId qpid, PathId spid, 
                           PathToAlnMap &path2aln_map,
						   int direction, 
								  Param &param );
IntPair findOverlap( PathId qpid,
                     PathId spid,
                     PathToAlnMap &path2aln_map,
					 int direction, 
							Param &param );

KmerId getCenterKmer( PathId query_path, 
                      PathId sbjct_path, 
                      int overlap,
                      int direction, 
                      PathToAlnMap &path2aln_map,
							 Param &param );

bool __checkLatchableRead( int &start,
                                 std::vector<KmerId> &kmers, 
                                 std::vector<int> &pov_vec, 
                                 PathId &path, 
                                 PathToAlnMap &path2aln_map, 
                                 Param &param, 
								  int direction);

void getLatchReads( BitString *bstrs,
                    PathId *used_reads,
                    PathId query_path,
                    PairedPath &pair_path,
                    KmerToPathMap &pathid_map,
                    PathToAlnMap &path2aln_map,
                    PathIdSet &merged_paths,
                    InvertedIndex &iindex,
						   size_t latch_offset, Param &param);

/*  void getLatchReads( ReadIdArray &latch_reads, */
/* 					 std::vector<int> &read_inits, */
/* 					 BitString *bstrs, */
/* 					 PathId *used_reads, */
/* 					 PathId query_path,  */
/* 					 //PairedPath &pair_path, */
/* 					 PathId sbjct_path, */
/* 					 KmerToPathMap &pathid_map, */
/* 					 PathToAlnMap &path2aln_map,  */
/* 					 PathIdSet &merged_paths, */
/* 					 InvertedIndex &iindex, */
/* 					 size_t latch_offset, Param &param); */

std::list<PairedPath> makePairedPaths( PathIdArray &pids,
                                       PathId &query_path,
                                       ReadId *pairs,
                                       PathId *used_reads,
                                       char *strands,
                                       PathToAlnMap &path2aln_map,
											  Param &param );

 void orderPathsByOverlap( std::multimap<int, PairedPath> &orders, std::list<PairedPath> &ppaths, Param &param  );
 void orderPathsBySupport( std::multimap<int, PairedPath> &orders, std::list<PairedPath> &ppaths, Param &param  );
 PairedPath determineBestPair( std::list<PairedPath> &ppaths, Param &param  );

std::vector<int> extractInits( ReadIdArray &rids, 
									  SpaPath *spath );

void stitchPair( PathId &query_path, 
                 PairedPath &pair_info,
                 KmerToPathMap &pathid_map,
                 PathToAlnMap &path2aln_map,
                 BitString *bstrs,
                 char *strands,
                 ReadId *pairs,
                 PathId *used_reads,
                 PathIdSet &merged_paths,
                 PathIdSet &success_paths,
                 InvertedIndex &iindex,
						Param &param );

bool linkPairedPath( PathId &query_path,
                     PairedPath &max_pair,
                     KmerToPathMap &pathid_map,
                     PathToAlnMap &path2aln_map, 
                     BitString *bstrs,
                     char *strands,
                     ReadId *pairs, 
                     PathId *used_reads,
                     int nreads, 
                     PathIdSet &merged_paths,
                     PathIdSet &success_paths,
                     InvertedIndex &iindex,
                     size_t latch_offset,
							Param &param );

void extendShortOverlapPaths( PathToAlnMap &path2aln_map, 
                       BitString *bstrs,
					   char *strands,
                       ReadId *pairs, 
                       PathId *used_reads,
					   int nreads, 
							  InvertedIndex &iindex,
									 Param &param );

 ScoreSummary getPathDepth(MSA &msa, Param &param);
 
 void printPositions( std::vector<int> &pos_vec );

ReadIdArray dropSelfPairs( ReadIdArray &orids,
                           PathId pid,
                           PathId *used_reads,
                           ReadId *pairs );


 void writeConsensus( std::fstream &out, SpaPath *spath , int count );
 void writePlacement( std::fstream &out, SpaPath *spath, int count );
 void writeStatistic( std::fstream &out, SpaPath *spath, int count );
 void writeAlignment( std::fstream &out, SpaPath *spath, int count, BitString *bstrs );
 void writeProfile( std::fstream &out, SpaPath *spath, int count );
};

#endif /* ASSEMBLY_H_ */

