#include "preprocessor.h"

using namespace std;

int main( int argc, char** argv )
{
    getopts( argc, argv );    
    proceed();
    return 0;
}

void proceed()
{
    BwaReader bwa( input_file, pairend );
    bwa.setDirectory(outdir);
    bwa.parse();

    GenomeMatrix gmat( genome_file.c_str() );
    gmat.setDirectory(outdir);
    gmat.setBwaReader(&bwa);
    gmat.generate();
}

void getopts( int argc, char **argv )
{
    int c;
    opterr = 0;
	while ((c = getopt (argc, argv, "i:g:o:pvh")) != -1) {
        switch (c) {
        case 'i' : input_file  = optarg; break;
        case 'g' : genome_file = optarg; break;
        case 'o' : outdir      = optarg; break;
        case 'v' : verbose     = true; break;
		case 'p' : pairend     = true; break;
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

    if ( input_file == "" ) {
        cerr << "BWA file is not given\n"; usage();
    }
    if ( genome_file == "" ) {
        cerr << "Genome file is not given\n"; usage();
    }
}

void usage()
{
    cerr << "Usage:program [options]\n"
         << "----------------------------------------------------------------------\n"
         << "  -i : string     required    BWA file\n"
		 << "  -g : string     required    genome information file\n"
         << "  -o : string     optional    output directory (default:.)\n"
		 << "  -p : flag       optional    pairend end reads\n"
         << "  -v : flag       optional    generate log message (default:off)\n"
         << "----------------------------------------------------------------------\n"
         << "e.g. bwa_handler -i file.bam -g genome_desc\n";
        
    exit(1);
}
