/*
  Youngik Yang
  Fri 2015-06-19 09:07:13 PM
 */

#ifndef _GEM_HANDLER_H__
#define _GEM_HANDLER_H__

#include "matrix.h"
#include "gem.h"
#include "mg.h"
#include "log.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <unistd.h>
#include <vector>
#include <typeinfo>

const double NOISE = 1.0e-12;
const int BOOT_ITER = 100;
const int DATA_RATE = 30;

/* Some global variable */
char *input = NULL;
char *lfile = NULL;
char *qfile = NULL;
std::string outdir = ".";
char   delim;
int nrow=0, ncol=0, ncpu=1, max_iter=MAX_ITER, dtype=0;
int    boot_iter = BOOT_ITER, sub_iter = BOOT_ITER;
bool   verbose = false, even = true;
bool   bootstrap = false, subsample=false;
double tolerance=TOLERANCE, pseudo=0, ts;
double percentage = DATA_RATE;


void parseOptions( int argc, char **argv );
void printCommand(int ac, char **av, std::ostream &out);
void printOptions();
void usage();
void loadMatrix( const char *input, Matrix &m, size_t nrow, size_t ncol, int dtype);
void report( GEM &g, MG &a, std::vector<std::string> &n);
void generateConfidenceInterval( Matrix &m, MG &g );

#endif
