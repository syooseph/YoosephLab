/*==============================================================================
 * Preprocessor
 * Author: Youngik Yang
 * Youngik Yang 
 * Created: Thu 2015-08-20 01:03:49 PM
 * Last modified: 
 *------------------------------------------------------------------------------
 * Input: BWA file in SAM format, genome description file
 * Output: Matrix (row: reads, column: genomes), preliminary abundance results
 *==============================================================================*/

#ifndef __PREPROCESSOR_H__
#define __PREPROCESSOR_H__

#include "genome_matrix.h"

/* global variables */
std::string input_file;
std::string genome_file;
std::string outdir = ".";
bool gzip_file     = false;
bool verbose       = false;
bool pairend       = false;

/* function prototypes */
void proceed();	 
void getopts( int argc, char **argv );
void usage();

#endif
