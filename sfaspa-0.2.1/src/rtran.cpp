/** 
 * @file rtran.cpp
 * @date Thu 2012-01-12 01:14:50 PM
 * @author Youngik Yang
 * @version 0.001
 * @brief Main function for reverse translation
 */

#include "rtran.h"


int main( int argc, char* argv[] ) 
{   
	std::string faa_file;
	std::string ffn_file;
	std::string place_file;
	std::string out_dir;
    int ncpus = 1;
    int gap_ext = GAP_EXT;
    int gap_open = GAP_OPEN;
    // int lower_band = LOWER_BAND;
    // int upper_band = UPPER_BAND;
    double band_ratio = BAND_RATIO;
    bool banded    = false;
    bool align     = false;
    bool profile   = false;
	bool verbose   = false;

    args::parse_cmd_args_postspa( argc, 
                                  argv, 
                                  faa_file,
                                  ffn_file,
                                  place_file,
                                  out_dir,
                                  ncpus,
                                  gap_ext,
                                  gap_open,
                                  //lower_band,
                                  //upper_band,
                                  band_ratio,
                                  banded,
                                  align,
                                  profile,
                                  verbose );

    args::printCommand(argc, argv, std::cout); 

    Reverser reverser( ffn_file, faa_file, place_file, out_dir, ncpus, align, profile, verbose );
    
    reverser.setGapPenalties( gap_ext, gap_open );

    if ( banded ) {
        reverser.setBandedAlignment(true);
        //reverser.setAlignmentBands( lower_band, upper_band );
        reverser.setBandRatio( band_ratio );
    }

    reverser.run();

    return 0;
}
