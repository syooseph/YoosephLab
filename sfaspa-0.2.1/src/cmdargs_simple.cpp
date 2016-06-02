#include "cmdargs.h"

void args::checkKmerSize(int k)
{
#if LONGKMER == 1
    if ( k > 12 ) {
        std::cerr << "[Error] Maximum k-mer size is 12\n";
        exit(EXIT_FAILURE);
    } 
#else
    if ( k > 6 ) {
        std::cerr << "[Error] Maximum k-mer size of the current binary is 6\n";
        std::cerr << "For the bigger kmer size, recompile with longer kmer support\n";
        exit(EXIT_FAILURE);
    }
#endif
}

void args::loadFileList(const char *listfile, std::vector<std::string> &input_files) 
{
    std::ifstream inf(listfile, std::ifstream::in );
    
    if ( !inf ) {
        std::cerr << "[Error] Can't open " << listfile << "\n";
        exit(EXIT_FAILURE);
    }
    
    std::string file;
    
    inf >> file;	
    while ( !inf.eof() ) {
        input_files.push_back( file );
        inf >> file;
    }
    inf.close();
}

void args::printCommand(int ac, char **av, std::ostream &out) 
{
    out << "Command: ";
    for ( int i = 0; i < ac; i++ )
        out << av[i] << " ";
    out << "\n\n";
}

void args::parse_cmd_args_postfgs(int ac, char **av, 
                            std::vector<std::string> &input_files, 
                            std::string &fgs_out,
                            std::string &fgs_ffn,
                            std::string &fgs_faa,
                            std::string &out_dir,
                            bool &verbose)
{
    try {
        std::string file_list;
        
        po::options_description generic("Generic options");
        generic.add_options()
            ("help,h", "produce help message")
            ;
        
        po::options_description list("File options");
        list.add_options()
            ("out-dir,t", po::value<std::string>(&out_dir)->default_value(OUTDIR), "output directory")
            //("input-list,l", po::value<std::string>(&file_list), "input file list")
            ("input-file,i", po::value< std::vector<std::string> >(), "input file")
            ("out-file,o", po::value<std::string>(&fgs_out), "FragGeneScan out file")
            ("ffn-file,n", po::value<std::string>(&fgs_ffn), "FragGeneScan ffn file")
            ("faa-file,a", po::value<std::string>(&fgs_faa), "FragGeneScan faa file")
            ("verbose,V", "verbose switch on")
            ;
		
        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden("Other options");
        hidden.add_options()
            ("input-file(s)", po::value< std::vector<std::string> >(), "input file")
            ;
		
        po::options_description cmdline_options;
        cmdline_options.add(generic).add(list).add(hidden);

        po::options_description visible("Options");
        visible.add(generic).add(list);
        
        po::positional_options_description p;
        p.add("input-file", -1);
        
        po::variables_map vm;
        store(po::command_line_parser(ac, av).
              options(cmdline_options).positional(p).run(), vm);
        notify(vm);


        if (vm.count("verbose")) 
            verbose = true;

        if (vm.count("help")) {
            std::cout << "\n" << USAGE_POSTFGS << "\n\n";
            std::cout << visible << "\n";
            exit(EXIT_FAILURE);
        }

        if ( ! vm.count("out-file") ||
             ! vm.count("ffn-file") ||
             ! vm.count("faa-file") ) {
            std::cerr << "[Error] FGS files must be provided\n";
            std::cout << "\n" << USAGE_POSTFGS << "\n\n";
            std::cout << visible;
            exit(EXIT_FAILURE);
        }

        if ( ! vm.count("input-list") && 
             ! vm.count("input-file") ) {
            std::cerr << "[Error] Input reads must be given\n";
            std::cout << "\n" << USAGE_POSTFGS << "\n\n";
            std::cout << visible;
            exit(EXIT_FAILURE);
        }
			
        if ( vm.count("input-list")) {
            file_list = vm["input-list"].as<std::string>();
            loadFileList(file_list.c_str(), input_files);
        }
        if (vm.count("input-file")) {
            input_files = vm["input-file"].as< std::vector<std::string> >();
        }

    }
    catch(std::exception& e) {
        std::cerr << e.what() << "\n";
        exit(EXIT_FAILURE);
    }    

}

void args::parse_cmd_args_parts(int ac, 
                                char **av, 
                                std::vector<std::string> &input_files )
{
    try {
        std::string file_list;

        /* Declare a group of options that will be allowed only on command line */
        //po::options_description generic("Generic options");
        po::options_description generic("");
        generic.add_options()
            ("help,h", "This help message")
            ("input-file,i", po::value< std::vector<std::string> >(), "Input file(s)")
            ("input-list,l", po::value<std::string>(&file_list), "File of multiple sequence locations")
            ;
			
        /* Hidden options, will be allowed both on command line and
         * in config file, but will not be shown to the user. */
        po::options_description hidden("Other options");
        hidden.add_options()
            ("input-file(s)", po::value< std::vector<std::string> >(), "input file")
            ;
		
        po::options_description cmdline_options;
        cmdline_options.add(generic).add(hidden);

        po::options_description visible("Options");
        visible.add(generic);
        
        po::positional_options_description p;
        p.add("input-file", -1);
        
        po::variables_map vm;
        store(po::command_line_parser(ac, av).
              options(cmdline_options).positional(p).run(), vm);
        notify(vm);

        if (vm.count("help")) {
            std::cout << "\n" << USAGE_PART << "\n\n";
            std::cout << visible << "\n";
            exit(EXIT_FAILURE);
        }

        if (vm.count("version")) {
            std::cout << VERSION << "\n";
            exit(EXIT_FAILURE);
        }

        if ( ! vm.count("input-list") && 
             ! vm.count("input-file") ) {
            std::cerr << "[Error] Input file(s) must be given\n";
            std::cout << "\n" << USAGE_PART << "\n\n";
            std::cout << visible;
            exit(EXIT_FAILURE);
        }
			
        if ( vm.count("input-list")) {
            file_list = vm["input-list"].as<std::string>();
            loadFileList(file_list.c_str(), input_files);
        }
        if (vm.count("input-file")) {
            input_files = vm["input-file"].as< std::vector<std::string> >();
        }
    }
    catch(std::exception& e) {
        std::cerr << e.what() << "\n";
        exit(EXIT_FAILURE);
    }    
}

void args::parse_cmd_args_prepare(int ac, 
                            char **av, 
                            int &k, 
                            std::vector<std::string> &input_files, 
                            std::string &outdir,
                            //bool &pair_flag,
                            //int &ncpu,
                            int &npart,
                                  bool &reverse_sfa,
                            bool &sfa_only,
                            bool &verbose)
{
    try {
        std::string file_list;

        /* Declare a group of options that will be allowed only on command line */
        //po::options_description generic("Generic options");
        po::options_description generic("");
        generic.add_options()
            //("version,v", "print version string")
            ("help,h", "This help message")
            ("kmer,k", po::value<int>(&k)->default_value(KMER_SIZE), "Size of k-mer")
            ("input-file,i", po::value< std::vector<std::string> >(), "Input file(s)")
            ("input-list,l", po::value<std::string>(&file_list), "File of multiple sequence locations")
            ("out-dir,o", po::value<std::string>(&outdir)->default_value(OUTDIR), "Output directory")
            ("parts,s", po::value<int>(&npart)->default_value(1), "no. of suffix array partitions")
            ;
			
        /* Hidden options, will be allowed both on command line and
         * in config file, but will not be shown to the user. */
        po::options_description hidden("Other options");
        hidden.add_options()
            ("full-help,H", "Full program options description")
            ("verbose,V", "verbose switch on")
            ("reverse-sfa,R", "Generate reverse suffix array")
            //("input-file(s)", po::value< std::vector<std::string> >(), "input file")
            ("suffix-only,S", "Do not generate graph & index")
            ;
		
        po::options_description cmdline_options;
        //cmdline_options.add(generic).add(parameter).add(list).add(hidden);
        cmdline_options.add(generic).add(hidden);

        po::options_description visible("Options");
        //visible.add(generic).add(parameter).add(list);
        visible.add(generic);

        po::options_description fullopts("Options");
        fullopts.add(generic).add(hidden);
        
        po::positional_options_description p;
        p.add("input-file", -1);
        
        po::variables_map vm;
        store(po::command_line_parser(ac, av).
              options(cmdline_options).positional(p).run(), vm);
        notify(vm);

        if (vm.count("help")) {
            std::cout << "\n" << USAGE_PREPARE << "\n\n";
            std::cout << visible << "\n";
            exit(EXIT_FAILURE);
        }

        if (vm.count("full-help")) {
            std::cout << "\n" << USAGE_ASSEM << "\n\n";
            std::cout << fullopts << "\n";
            exit(EXIT_FAILURE);
        }

        if (vm.count("version")) {
            std::cout << VERSION << "\n";
            exit(EXIT_FAILURE);
        }

        if ( ! vm.count("input-list") && 
             ! vm.count("input-file") ) {
            std::cerr << "[Error] Input file(s) must be given\n";
            std::cout << "\n" << USAGE_PREPARE << "\n\n";
            std::cout << visible;
            exit(EXIT_FAILURE);
        }
			
        if ( vm.count("input-list")) {
            file_list = vm["input-list"].as<std::string>();
            loadFileList(file_list.c_str(), input_files);
        }
        if (vm.count("input-file")) {
            input_files = vm["input-file"].as< std::vector<std::string> >();
        }

        if (vm.count("suffix-only"))
            sfa_only = true;

        if (vm.count("verbose"))
            verbose = true;

        if (vm.count("reverse-sfa"))
            reverse_sfa = true;

        /* Error check */
        checkKmerSize(k);
    }
    catch(std::exception& e) {
        std::cerr << e.what() << "\n";
        exit(EXIT_FAILURE);
    }    
}

	
void args::parse_cmd_args_assembly( int ac,
                                    char **av,
                                    Param &param)
{
    //std::string merge_score_str = std::to_string( MERGE_SCORE );
    char temp[256];
    std::sprintf(temp, "%.1f", MERGE_SCORE);
    std::string merge_score_str = temp;
    std::sprintf(temp, "%.1f", MERGE_SHORT_RATIO);
    std::string merge_short_ratio_str = temp;
    std::sprintf(temp, "%.1f", EXTEND_SCORE);
    std::string extend_score_str = temp;
    std::sprintf(temp, "%.1f", READ_ALIGN_SCORE);
    std::string read_align_score_str = temp;
    std::sprintf(temp, "%.1f", READ_ALIGN_RATIO);
    std::string read_align_ratio_str = temp;

    try {
        std::string file_list;

        po::options_description generic;
        generic.add_options()
            ("help,h", "This help message")
            ("pair-end,M", "Paired-end reads [Input reads needs to be Interlaced]")
            ("input-file,i", po::value< std::vector<std::string> >(), "Sequence file(s)")
            //("input-list,l", po::value<std::string>(&file_list), "File of multiple sequence locations")
            ("output-dir,o", po::value<std::string>(&param.out_dir)->default_value(OUTDIR), "Output directory")
            ("kmer-size,k", po::value<int>(&param.kmer_size )->default_value(KMER_SIZE), "Size of k-mer (node)")
            ("nparts,s", po::value<int>(&param.nparts )->default_value(NUM_PARTS), "No. of suffix array partitions")
            ("ncpus,n", po::value<int>(&param.ncpus )->default_value(NUM_CPUS), "No. of CPUs")
            ("seed-coverage,S", po::value<int>(&param.min_seed)->default_value(MIN_SEED), "Minimum seed coverage")
            ("read-support,r", po::value<int>(&param.min_share)->default_value(SHARED_READS), "Minimum read support")
            ("overlap-length,b", po::value<int>(&param.back_trace)->default_value(BACK_TRACE), "Minimum overlap length for suffix array search")
            ("seed-reuse,u", po::value<bool>(&param.seed_reuse)->default_value(SEED_REUSE), "Reuse seed k-mers found in previous paths for seeding paths")
            ("profile,P", "Write a profile file")
            ("alignment,A", "Write an alignment file")
            ("version,v", "Print SPA version")
            ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden("Advanced options");
        hidden.add_options()
            ("full-help,H", "Full program options description")
            //("input-file(s)", po::value< std::vector<std::string> >(), "input file")
            ("verbose,V", po::value<int>(&param.verbose)->default_value(VERBOSITY), "Verbose option (0:silent, 1:brief, 2:detailed, 3:very wordy)")
            //("summary,B", po::value<bool>(&param.summary_flag)->default_value(false), "Print summary and information during assembly")
            ("output-all,a", po::value<bool>(&param.output_all)->default_value(OUTPUT_ALL), "Generate outputs of each assembly stage")
            ("gap-open,g", po::value<int>(&param.gap_open)->default_value(GAP_OPEN), "Gap open penalty")
            ("gap-ext,x", po::value<int>(&param.gap_ext)->default_value(GAP_EXT), "Gap extension penalty")
            //("trim-graph,T", po::value<bool>(&param.trim_flag)->default_value(false), "Trim low coverage vertices/edges")
            //("min-coverage,c", po::value<int>(&param.min_depth)->default_value(MIN_DEPTH), "Minimum kmer coverage for trimming graph")

            //("path-length,l", po::value<int>(&param.min_length)->default_value(MIN_LENGTH), "Minimum path length in path extraction stage")
            //("dump-flag,D", po::value<bool>(&param.dump_flag)->default_value(false), "Dump output as binary")
            //("extract-paths,E", po::value<bool>(&param.path_flag)->default_value(true), "Initial path search")
            //("cluster-paths,C", po::value<bool>(&param.merge_flag)->default_value(true), "Merge similar paths")
            //("recluster-paths,G", po::value<bool>(&param.recluster_flag)->default_value(false), "Clustering similar paths again after path extensions")
            //("extend-paths,J", po::value<bool>(&param.extend_flag)->default_value(true), "Connect paths")
            //("place-reads,L", po::value<bool>(&param.place_flag)->default_value(true), "Place reads to paths")
            //("recruit-reads,R", po::value<bool>(&param.recruit_flag)->default_value(true), "Recruit unassigned reads")
            //("report-outputs,O", po::value<bool>(&param.report_flag)->default_value(true), "Make SPA output")

            //("cluster-expansion", po::value<bool>(&param.extend_cluster)->default_value(true), "Allow pivot sequences to be elongated")            
            //("align-base-first", po::value<bool>(&param.align_base_first)->default_value(true), "Perform simple alignment by base comparsion before full alignment")            
            //("banded-alignment", po::value<bool>(&param.banded_align)->default_value(BAND_ALIGN), "Perform banded alignment")
            //("band-ratio", po::value<double>(&param.band_ratio)->default_value(BAND_RATIO, boost::lexical_cast<std::string>(boost::format("%.2f") % BAND_RATIO)), "Max band ratio w.r.t. a sequence length")
            ("insert-size", po::value<int>(&param.insert_size)->default_value(INSERT_SIZE), "Sequence library inner insert size")
            ("insert-sd", po::value<double>(&param.insert_sd)->default_value(INSERT_SD), "Standard deviation of insert size")
            //("merge-filter-kmer", po::value<int>(&param.merge_filter_kmer)->default_value(MERGE_FILTER_KMER), "Initital k-mer size for clustering")
            //("merge-shared-nkmer", po::value<int>(&param.merge_shared_nkmer)->default_value(MERGE_SHARED_NKMER), "Minimum no. of shared kmers for k-mer filtering")
            //("merge-expand-nbase", po::value<int>(&param.merge_expand_nbase)->default_value(MERGE_EXPAND_NBASE), "No. of extra bases to make filter region wider")
            //("merge-filter-score", po::value<double>(&param.merge_filter_score)->default_value(MERGE_FILTER_SCORE, boost::lexical_cast<std::string>(boost::format("%.2f") % MERGE_FILTER_SCORE)), "Filtering score for clustering")
            //("merge-score", po::value<double>(&param.merge_score)->default_value(MERGE_SCORE, boost::lexical_cast<std::string>(boost::format("%.2f") % MERGE_SCORE)), "Alignment score for clustering")
            ("merge-score", po::value<double>(&param.merge_score)->default_value(MERGE_SCORE, merge_score_str), "Alignment score for clustering")
            //("merge-fixed-kmer", po::value<bool>(&param.merge_fixed_kmer)->default_value(false), "Use fixed size of filter k-mer for similar path search during clustering")
            //("merge-short-ratio", po::value<double>(&param.merge_short_ratio)->default_value(MERGE_SHORT_RATIO, boost::lexical_cast<std::string>(boost::format("%.2f") % MERGE_SHORT_RATIO)), "Minium aligned portion of shorter sequence")
            ("merge-short-ratio", po::value<double>(&param.merge_short_ratio)->default_value(MERGE_SHORT_RATIO, merge_short_ratio_str), "Minium aligned portion of shorter sequence")
            //("extend-anchor-kmer", po::value<int>(&param.extend_anchor_kmer)->default_value(EXTEND_ANCHOR_KMER), "anchor Kmer size of path extension")
            //("extend-anchor-mink", po::value<int>(&param.extend_anchor_mink)->default_value(EXTEND_ANCHOR_MINK), "minimum shared anchor kmers")
            //("extend-filter-kmer", po::value<int>(&param.extend_filter_kmer)->default_value(EXTEND_FILTER_KMER), "Kmer size for overlap extension")
            //("extend-filter-score", po::value<double>(&param.extend_filter_score)->default_value(EXTEND_FILTER_SCORE, boost::lexical_cast<std::string>(boost::format("%.2f") % EXTEND_FILTER_SCORE)), "Filtering score for overlap extension")
            //("extend-off-nbase", po::value<int>(&param.extend_off_nbase)->default_value(EXTEND_OFF_NBASE), "Off diagonal allowance for path extension region")
            //("extend-score", po::value<double>(&param.extend_score)->default_value(EXTEND_SCORE, boost::lexical_cast<std::string>(boost::format("%.2f") % EXTEND_SCORE)), "Alignment score for overlap extension")
            ("extend-score", po::value<double>(&param.extend_score)->default_value(EXTEND_SCORE, extend_score_str), "Alignment score for overlap extension")
            //("long-overlap-only", po::value<bool>(&param.long_overlap_only)->default_value(false), "Long overlap extension only")
            //("short-overlap-only", po::value<bool>(&param.short_overlap_only)->default_value(false), "Short overlap extension only")
            //("read-overlap-only", po::value<bool>(&param.read_overlap_only)->default_value(false), "Read overlap extension only")
            //("read-bridge-extend", po::value<bool>(&param.read_bridge_extend)->default_value(true), "Path extension with bridging reads")
            //("recruit-filter-kmer", po::value<int>(&param.recruit_filter_kmer)->default_value(RECRUIT_FILTER_KMER), "Kmer size for read recruitment")
            //("recruit-min-filter", po::value<int>(&param.recruit_min_filter)->default_value(RECRUIT_MIN_FILTER), "Minimum kmer match for read recruitment")
            //("recruit-filter-score", po::value<double>(&param.recruit_filter_score)->default_value(RECRUIT_FILTER_SCORE, boost::lexical_cast<std::string>(boost::format("%.2f") % RECRUIT_FILTER_SCORE)), "Filtering score for read recruitment")
            //("recruit-score", po::value<double>(&param.read_align_score)->default_value(READ_ALIGN_SCORE, boost::lexical_cast<std::string>(boost::format("%.2f") % READ_ALIGN_SCORE)), "Alignment score for read recruitment")
            ("recruit-score", po::value<double>(&param.read_align_score)->default_value(READ_ALIGN_SCORE, read_align_score_str), "Alignment score for read recruitment")
            ("recruit-length", po::value<double>(&param.read_align_ratio)->default_value(READ_ALIGN_RATIO, read_align_ratio_str), "Alignment length ratio for read recruitment")
            //("recruit-by-align", po::value<bool>(&param.recruit_by_align)->default_value(false), "Recruit reads by alignment instead of anchoring")
            //("use-identity", po::value<bool>(&param.identity_flag)->default_value(false), "Use indentity instead of positive score")
            //("use-percentile", po::value<bool>(&param.percent_flag)->default_value(false), "Use percentile during seed selection, instead of fixed coverage value")
            //("percentile", po::value<double>(&param.kmer_percentile)->default_value(KMER_PERCENTILE, boost::lexical_cast<std::string>(boost::format("%.2f") % KMER_PERCENTILE)), "Top x% percentile of seed kmers")
            //("par-search", po::value<bool>(&param.par_search)->default_value(true), "Best neighboring node search in parallel")
            //("skip-fail", po::value<bool>(&param.skip_fail)->default_value(true), "Avoid seed kmers from failed paths")
            //("try-later", po::value<bool>(&param.try_later)->default_value(true), "Try failed seed later if skipped")
            //("use-all-fail", po::value<bool>(&param.use_all_fail)->default_value(true), "Try all filed nodes as seeds")
            //("strict-length", po::value<bool>(&param.strict_length)->default_value(false), "Path length must be at least min. length")
            //("trimmed-short", po::value<bool>(&param.short_trim_len)->default_value(true), "Allow shortened path after trimming?")
            //("node-check", po::value<bool>(&param.path_node_check)->default_value(true), "Check existence of zero coverage path node")
            ("extend-first", po::value<bool>(&param.extend_first)->default_value(true), "Path extension before clustering")
            //("check-stop", po::value<bool>(&param.check_stop)->default_value(false), "Check stop condition conflict during extension")
            //("drop-dup-path", po::value<bool>(&param.drop_dup_path)->default_value(true), "Remove subpath duplicates")
            //("check-cycle-entry", po::value<bool>(&param.check_cycle_entry)->default_value(false), "Check cycle entry points")
            ("extend-length", po::value<int>(&param.extend_length)->default_value(EXTEND_LENGTH), "Mininum length of long overlapping extension")
            ("suffix-overlap", po::value<int>(&param.suffix_overlap)->default_value(SUFFIX_OVERLAP), "Mininum length of short overlapping extension by suffix array support")
            ("bridge-overlap", po::value<int>(&param.bridge_overlap)->default_value(BRIDGE_OVERLAP), "Mininum overlapping length between a read and end path to find bridging reads")
            ("min-bridge-reads", po::value<int>(&param.min_bridge_reads)->default_value(MIN_BRIDGE_READS), "Mininum number of bridging reads")
            //("ignore-stop", po::value<bool>(&param.ignore_stop)->default_value(false), "Ignore stop codon while path extension")
            //("replace-stop", po::value<bool>(&param.replace_stop)->default_value(true), "Replace stop codon with other amino acid during consensus call, if possible")
            //("clip-stop", po::value<bool>(&param.clip_stop)->default_value(true), "Clip sequence away after stop codon")
            //("release-path", po::value<bool>(&param.release_path)->default_value(true), "Release path after initial path search to avoid excessive memory usage")
            ("ljoin-pe-support", po::value<bool>(&param.ljoin_pe_support)->default_value(LONG_OVERLAP_PATH_PAIREND_SUPPORT), "Should have pair-end support for long overlapping path extension?")
            ("sjoin-pe-support", po::value<bool>(&param.sjoin_pe_support)->default_value(SHORT_OVERLAP_PATH_PAIREND_SUPPORT), "Should have pair-end support for short overlapping path extension?")
            ("rjoin-pe-support", po::value<bool>(&param.rjoin_pe_support)->default_value(READ_BRIDGE_PATH_PAIREND_SUPPORT), "Should have pair-end support for read bridging path extension?")
            ("min-pair-reads", po::value<int>(&param.min_pair_reads)->default_value(MIN_PAIR_READS), "Minimum number of pairend reads")
            //("debug-id", po::value<int>(&param.debug_id), "Object ID for debugging purpose only")
            //("debug-flag", po::value<bool>(&param.debug_flag), "Print out debugging portions")
            //("exact-match-only", po::value<bool>(&param.exact_match_only)->default_value(true), "Exact matching reads are placed during path extraction")
            ;


        po::options_description cmdline_options;
        cmdline_options.add(generic).add(hidden);

        po::options_description visible("Options");
        visible.add(generic);

        po::options_description fullopts("Options");
        fullopts.add(generic).add(hidden);
                                         
        po::positional_options_description p;
        p.add("input-file", -1);

        po::variables_map vm;
        store(po::command_line_parser(ac, av).
              options(cmdline_options).positional(p).run(), vm);
        notify(vm);


        if (vm.count("version")) {
            std::cout << VERSION << "\n";
            exit(EXIT_FAILURE);
        }
        if (vm.count("help")) {
            std::cout << "\n" << USAGE_ASSEM << "\n\n";
            std::cout << visible << "\n";
            exit(EXIT_FAILURE);
        }
        if (vm.count("full-help")) {
            std::cout << "\n" << USAGE_ASSEM << "\n\n";
            std::cout << fullopts << "\n";
            exit(EXIT_FAILURE);
        }

        if (vm.count("version")) {
#ifdef COMMIT
            printf("Version: %s\n", COMMIT);  // git commit version
#else
            //printf("Version: %s\n", VERSION); // preset variable
            std::cout << "Version: " << VERSION << "\n";            
#endif
            exit(EXIT_FAILURE);
        }

        if ( ! vm.count("input-list") && !vm.count("input-file") ) {
            std::cerr << "\n[Error] sequence file(s) must be given\n";
            std::cout << "\n" << USAGE_ASSEM << "\n\n";
            std::cout << visible;
            exit(EXIT_FAILURE);
        }

        //=========================================================
        // I am going to remove uniq_seed parameter later on
        // But I am going to keep now since it uses in main program
        //=========================================================
        //param.uniq_seed = !param.seed_reuse;

        if (vm.count("pair-end"))
            param.pair_flag = true;
        if (vm.count("profile"))
            param.profile_flag = true;
        if (vm.count("alignment"))
            param.align_flag = true;
        // if (vm.count("verbose"))
        //     param.verbose = true;

        // if (vm.count("super-verbose"))
        //     param.super_verbose = param.verbose = true;

        if ( vm.count("input-list")) {
            file_list = vm["input-list"].as<std::string>();
            loadFileList(file_list.c_str(), param.input_files);
        }
        if (vm.count("input-file")) {
            param.input_files = vm["input-file"].as< std::vector<std::string> >();
        }

        
        checkKmerSize(param.kmer_size);

        //if ( param.verbose == 1 ) param.summary_flag = 1;
    }
    catch(std::exception& e) {
        std::cerr << e.what() << "\n";
        exit(EXIT_FAILURE);
    }

}


void args::parse_cmd_args_postspa( int ac, 
                                   char **av, 
                                   std::string &faa_file,
                                   std::string &ffn_file,
                                   std::string &place_file,
                                   std::string &out_dir,
                                   int &ncpus,
                                   int &gap_ext,
                                   int &gap_open,
                                   //int &lower_band,
                                   // int &upper_band,
                                   double &band_ratio,
                                   bool &banded,
                                   bool &align,
                                   bool &profile,
                                   bool &verbose )
{

    std::string band_ratio_str = std::to_string(BAND_RATIO);
    try {

        po::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "Print version string")
            ("help,h", "Produce help message")
            ;

        po::options_description list("File options");
        list.add_options()
            ("faa-file,a", po::value<std::string>(&faa_file), "ffa file")
            ("ffn-file,d", po::value<std::string>(&ffn_file), "ffn file")
            ("place-file,p", po::value<std::string>(&place_file), "Placement file")
            ("out-dir,o", po::value<std::string>(&out_dir)->default_value("."), "Output directory")
            ("ncpus,n", po::value<int>(&ncpus)->default_value(1), "No of CPUs")
            ("align,A", "Write an alignment output file")
            ("profile,P", "Write a profile output file")
            ;
		
			
        /* Hidden options, will be allowed both on command line and
         * in config file, but will not be shown to the user. */
        po::options_description hidden("Other options");
        hidden.add_options()
            ("full-help,H", "Full program options description")
            ("verbose,V", "Verbose")
            ("banded-alignment,B", po::value<bool>(&banded)->default_value(BAND_ALIGN), "Perform banded alignment")
            ("gap-ext,x", po::value<int>(&gap_ext)->default_value(GAP_EXT), "Gap extension penalty")
            ("gap-open,g", po::value<int>(&gap_open)->default_value(GAP_OPEN), "Gap open penalty")
            // ("lower-band,L", po::value<int>(&lower_band)->default_value(LOWER_BAND), "Lower band limit")
            // ("upper-band,U", po::value<int>(&upper_band)->default_value(UPPER_BAND), "Upper band limit")
            //("band-ratio,R", po::value<double>(&band_ratio)->default_value(BAND_RATIO, boost::lexical_cast<std::string>(boost::format("%.2f") % BAND_RATIO)), "Max band ratio w.r.t. a sequence length")
            ("band-ratio,R", po::value<double>(&band_ratio)->default_value(BAND_RATIO, band_ratio_str), "Max band ratio w.r.t. a sequence length")
            ;

        po::options_description cmdline_options;
        cmdline_options.add(generic).add(list).add(hidden);

        po::options_description visible("Options");
        visible.add(generic).add(list);

        po::options_description fullopts("Options");
        fullopts.add(generic).add(list).add(hidden);

        po::variables_map vm;
        store(po::command_line_parser(ac, av).
              options(cmdline_options).run(), vm);
        notify(vm);

        if (vm.count("verbose"))
            verbose = true;
        // if (vm.count("no-banded"))
        //     banded = false;
        if (vm.count("align"))
            align = true;
        if (vm.count("profile"))
            profile = true;

        if (vm.count("help")) {
            std::cout << "\n" << USAGE_POST_SPA << "\n\n";
            std::cout << visible << "\n";
            exit(EXIT_FAILURE);
        }

        if (vm.count("full-help")) {
            std::cout << "\n" << USAGE_POST_SPA << "\n\n";
            std::cout << fullopts << "\n";
            exit(EXIT_FAILURE);
        }

        if (vm.count("version")) {
            std::cout << VERSION << "\n";
            exit(EXIT_FAILURE);
        }

        if ( ! vm.count("faa-file") ||
             ! vm.count("ffn-file") ||
             ! vm.count("place-file") ) { 
            std::cerr << "[Error] .ffa .ffn .place files must be given\n";
            std::cout << "\n" << USAGE_POST_SPA << "\n\n";
            std::cout << visible;
            exit(EXIT_FAILURE);
        }
    }
    catch(std::exception& e) {
        std::cerr << e.what() << "\n";
        exit(EXIT_FAILURE);
    }    
}

