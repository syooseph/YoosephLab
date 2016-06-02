/*
 * Wed 2015-06-24 03:31:51 PM
 * Last modified: Mon 2015-08-03 01:39:33 PM
 * Youngik Yang
 */

#ifndef __METAGENOME_H__
#define __METAGENOME_H__

#include <stdint.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <cassert>

using SizeMap = std::unordered_map<std::string, uint32_t>;
using RateMap = std::unordered_map<std::string, double>;

class MG
{
 private:
	SizeMap size_map;
	RateMap abud_map;
	
	void loadGenomeSizes( const char *file );

 public:
	MG( const char *file);
	~MG();

	/* Compute genome relative abundances */
	void computeAbundances(double *coefficient, std::vector<std::string> &names);

	/* Getters */
	SizeMap getGenomeSizes() { return size_map; }
	RateMap getAbundances()  { return abud_map; } 
};

#endif
