#include "gem_handler.h"

using namespace std;

int main( int argc, char **argv )
{
	ts = mytime();
	
	//------------------
    // Disable buffering
    //------------------
    setbuf(stdout, NULL);
    

    parseOptions( argc, argv );

    //-------------
    // Input matrix
    //-------------
    Matrix m(input, delim);
	std::vector<string> colnames = m.getColNames();    
	string message = "Matrix loaded";
#if SPARSE == 1
	message = "Sparse " + message;
#endif
    log(message, ts);
	
	//-----------------
    // Matrix dimemsion
    //-----------------
    nrow = m.getRowSize();
    ncol = m.getColSize();

	//--------------------
    // Load genome lengths
    //--------------------
    MG mg(lfile);
    log( "Genome length loaded", ts );

    //-----------
    // Run Grammy
    //-----------
	GEM em(m.getMatrix(), nrow, ncol, ncpu, max_iter, tolerance, pseudo, even, verbose);
    em.run();
    log( "EM completed", ts );
	
    //-------------------------------
    // Get genome relative abundances
    //-------------------------------
    mg.computeAbundances( em.getCoefficient(), colnames );
    log( "Relative abundances computed", ts );
    
    //---------------------
    // Report Grammy result
    //---------------------
    report(em, mg, colnames);
    
    //------------------------------------------
    // Compute confidence interval by resampling
    //------------------------------------------
    if ( bootstrap || subsample ) {
        MG bmg = mg;
        generateConfidenceInterval( m, bmg);
    }

    log( "Done", ts);
    
    return 0;
}

void parseOptions( int argc, char **argv )
{
    int c;
    opterr = 0;
	while ((c = getopt (argc, argv, "hi:r:c:d:l:q:n:e:vmst:p:o:bB:P:")) != -1) {
        switch (c) {
        case 'i' : input     = optarg; break;
        case 'l' : lfile     = optarg; break;
        case 'q' : qfile     = optarg; break;
        case 'n' : ncpu      = atoi(optarg); break;
        case 'd' : dtype     = atoi(optarg); break;
        case 'e' : even      = (bool)atoi(optarg); break;
        case 'v' : verbose   = true; break;
        case 'm' : max_iter  = atoi(optarg); break;
        case 't' : tolerance = atof(optarg); break;
        case 'o' : outdir    = string(optarg); break;
        case 'b' : bootstrap = true; break;
        case 'B' : boot_iter = atoi(optarg); break;
        case 'P' : percentage  = atof(optarg); break;            
        case 'h' : usage(); 
        case '?':
            if (optopt == 'c')
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            else if (isprint (optopt))
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
                fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
            exit(1);
        default:
            usage();
        }
    }

    if ( input == NULL ) {
        cerr << "Input file is not given\n"; usage();
    }
    if ( lfile == NULL ) {
        cerr << "Genome length file is not given\n"; usage();
    }

    switch(dtype) {
    case 1 : delim = ' ';break;
    case 2 : delim = ','; break;
    default: delim = '\t';
    }

	printCommand(argc, argv, std::cout);
	printOptions();
}

void usage()
{
    cerr << "Usage:program [options]\n"
         << "----------------------------------------------------------------------\n"
         << "Input options:\n"
         << "  -i : string     required    input matrix file (csv format)\n"
         << "  -d : integer    optional    delimiter charcter of csv (0.tab, 1:space, 2:comma [default:0])\n"        
         << "  -l : string     required    genome length file\n"
         << "  -q : string     optional    read name file\n"        
         << "----------------------------------------------------------------------\n"        
         << "EM options:\n"
         << "  -n : integer    optional    number of cpu (default:1)\n"
         << "  -e : boolean    optional    initiaize mixing coefficients with equal weight (default:1)\n"
         << "  -m : integer    optional    maximum iteration (default:10000)\n"
         << "  -t : numeric    optional    convergency tolerance (default:1.0e-6)\n"
         << "----------------------------------------------------------------------\n"
         << "Other opitons:\n"
         << "  -v : flag       optional    verbose (default:off)\n"        
         << "  -o : string     optional    output directory (default:.)\n" 
         << "  -b : flag       optional    confidence interval by bootstrapping (default:off)\n"
         << "  -B : integer    optional    no. of bootstrapping (default:100)\n"
         << "  -P : numeric    optional    percentage of data to be sampled (default:30)\n";
    exit(1);
}

void printCommand(int ac, char **av, std::ostream &out) 
{
	std::cout << "Sparse matrix version:";
	SPARSE ? printf("Yes\n") : printf("No\n");
    char host[256];
    gethostname(host, 256);
    out << "Machine: " << host << "\n";
    out << "Command: ";
    for ( int i = 0; i < ac; i++ )
        out << av[i] << " ";
    out << "\n\n";
}

void printOptions()
{
	std::cout << "------------------------------\n";
	std::cout << "OPTIONS:\n";
	std::cout << "Convergence:" << std::scientific << tolerance << "\n";
	std::cout << "Max. iteration:" << max_iter << "\n";
	std::cout << "CPUs:" << ncpu << "\n";
	std::cout << "Boostrapping:";
	bootstrap ? printf("Yes\n") : printf("No\n");
	if ( bootstrap ) {
		std::cout << "\tIteration:" << boot_iter << "\n"
				  << "\tSample:" << percentage << "%\n";
	}
	std::cout << "------------------------------\n";
}

void report( GEM &em,
             MG  &mg,
             std::vector<string> &genomes)
{
    cout << "EM Convergency: " << em.converged() << "\n";
    cout << "EM Iterations: " << em.getIteration() << "\n";

	double *mixing = em.getCoefficient();
    string ofile = outdir + "/gabe_coefficients";
    ofstream out(ofile.c_str());
    if ( !out ) {
        cerr << "Can't open " << ofile << "\n"; exit (1);
    }
    
    for ( int i = 0; i < ncol; i++ ) 
        out << genomes[i] << "\t" << mixing[i] << "\n";
    out.close();
    cout << "EM mixing coefficients: " << ofile << "\n";

    RateMap abd_map = mg.getAbundances();
    ofile = outdir + "/gabe_abundances";
    out.open(ofile.c_str());
    if ( !out ) {
        cerr << "Can't open " << ofile << "\n"; exit (1);
    }
	
    for ( int i = 0; i < ncol; i++ ) 
        out << genomes[i] << "\t" << abd_map[genomes[i]] << "\n";
    out.close();
    cout << "Relative genome abundances: " << ofile << "\n";
    
    //-----------------
    // class assignment
    //-----------------
	std::vector<string> rnames(nrow, "");
	std::vector<size_t> edists(nrow, 0); // edit distances
    if ( qfile != NULL ) {
        ifstream in( qfile );
        if ( !in ) {
            cerr << "Can't open " << qfile << "\n"; exit (1);
        }
        size_t n = 0;
        string line;
        while(getline(in,line)) {
			if ( line[0] == '#' ) continue;
			std::istringstream iss(line);
			std::vector<string> c;
			std::string tok;
			while(iss>>tok) c.push_back(tok);
			assert(c.size()>=3);
            rnames[n] = c[1];
			edists[n] = stoul(c[2]);
            n++;
        }
        assert(n=nrow);
    }
	
    

	MatrixType *post = em.getPosterior();
    std::vector<int> MaxCol(nrow,0);
	std::vector<double> MaxVal(nrow,0);
	int i;
#pragma omp parallel for schedule(static, 1) if (ncpu>1) private(i) num_threads(ncpu)                    
    for ( i = 0; i < nrow; i++ ) {
        int max_pos = 0;
        double max_val = (*post)(i,0);
        for ( int j = 1; j < ncol; j++ ) {
			double v = (*post)(i,j);
            if ( v > max_val ) {
                max_pos = j;
                max_val = v;
            }
        }
		MaxCol[i] = max_pos;
		MaxVal[i] = max_val;
	}
	
    ofile = outdir + "/gabe_assignments";
    out.open(ofile.c_str());
    if ( !out ) {
        cerr << "Can't open " << ofile << "\n"; exit (1);
    }

	out << "Read_index\t";
	if ( qfile != NULL ) 
		out << "Read_name\tEdit_distance\t";
	out << "Genome_index\tGenome_name\tProbability\n";

    //-----------------
    // class assignment
    //-----------------
	for ( i = 0; i < nrow; i++ ) {
		int max_pos = MaxCol[i];
		out << i+1 << "\t";
		if ( qfile != NULL ) 
			out << rnames[i] << "\t" << edists[i] << "\t";
		out << max_pos+1 << "\t" << genomes[max_pos] << "\t" << MaxVal[i] << "\n";
	}
    out.close();
    cout << "Read assignments: " << ofile << "\n";
}

void generateConfidenceInterval( Matrix &m,
                                 MG &mg )
{
	std::vector<string> genomes = m.getColNames();

    if ( bootstrap )
        cout << "\nBootstrapping ...\n";
    else
        cout << "\nSubsampling ...\n";
        
    int niter = bootstrap ? boot_iter : sub_iter;
    double **abunds = new double *[niter];
    for ( int i = 0; i < niter; i++ )
        abunds[i] = new double[ncol];

    //--------------------
    // No. of sampled data
    //--------------------
    int nsample = nrow*percentage/100;    

	MatrixType sampled;
    //----------------------------------
    // Perform bootstrapping/subsampling
    //----------------------------------
    for ( int i = 0; i < niter; i++ ) {
		sampled = MatrixType(nsample, ncol);
        bootstrap ? m.resample(&sampled, nsample, i) : m.subsample(&sampled, nsample, i);
        GEM *gem = new GEM(&sampled, nsample, ncol, ncpu, max_iter, tolerance, pseudo, even, verbose);
        gem->run();

        bool con = gem->converged();
        size_t nit = gem->getIteration();

        mg.computeAbundances( gem->getCoefficient(), genomes );
        RateMap abd_map = mg.getAbundances();

        // Store genome abundances
        for ( int c = 0; c < ncol; c++ )
            abunds[i][c] = abd_map[genomes[c]];

        // Delete EM object
        delete gem;

        string message = bootstrap ? "#Bootstrap:" : "#Subsampling:";
        message += to_string(i+1);
        message += "\t(converged:" + to_string(con) + ", iteration:" + to_string(nit) + ")";
        log( message, ts );
    }
	
    string ofile = outdir + "/gabe_intervals";
    ofstream out(ofile.c_str());
    if ( !out ) {
        cerr << "Can't open " << ofile << "\n"; exit (1);
    }
    out << "# Genome abundance confidence interval (95%)\n";
    out << fixed; out.precision(8);
	std::vector<double> iters(niter,0.0);
    for ( int c = 0; c < ncol; c++ ) {
        for ( int i = 0; i < niter; i++ ) 
            iters[i] = abunds[i][c];
        sort( iters.begin(), iters.end() );
        size_t lq = round( iters.size() * 2.5/100.0 );
        size_t rq = round( iters.size() * 97.5/100.0 );

        double lb = lq > 0 ? iters[ lq-1 ] : iters[ 0 ];
        double rb = iters[ rq-1 ];
        out << genomes[c] << "\t[" << lb << "," << rb << "]\n";
    }
    out.close();

    cout << "Confidence interval: " << ofile << "\n";        
    
    for ( int i = 0; i < niter; i++ )
        delete[] abunds[i];
    delete[] abunds;
}
