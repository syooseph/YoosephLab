/** 
 * @file spamain.cpp
 * @date Thu Mon 2013-08-12 04:07:05 PM 
 * @date Modified: Thu 2014-02-27 11:41:13 AM
 * @author Youngik Yang
 * @version 0.002
 * @brief Main function of SPA assembler.
 * @details 
 * This is main handler of SPA assembler. 
 * It reads program options and run assembler.
 */

#include "spamain.h"

//====================================================================
// Main function of SPA
//====================================================================
int main( int argc, char* argv[] ) 
{	
    //-------------
    // No buffering
    //-------------
    setbuf(stdout, NULL); 

    //---------------------
    // Read program options
    //---------------------
    Param p;
    readParams(argc, argv, p);

    //--------------
    // Run assembler
    //--------------
    Assembler assem(p);
    assem.run();

    printElapsed(INIT_TIME, mytime(), "Done");
    printLocalTime();

    return 0;
}

//====================================================================
// Read program options from command line
//====================================================================
void readParams(int argc, char **argv, Param &p)
{
    //------------
    // Get options
    //------------
    args::parse_cmd_args_assembly( argc, argv, p );

    printLocalTime();

    //----------------------
    // Show current version
    //----------------------
// #ifdef COMMIT
//     printf("Version: %s\n", COMMIT);  // git commit version
// #else
    std::cout << "Version: " << VERSION << "\n";
#ifdef COMMIT
    printf("Git-tag: %s\n", COMMIT);  // git commit version
#endif

    //--------------
    // Print command
    //--------------
    args::printCommand(argc, argv, std::cout);
}

