/**
 * \file core.h 
 * \brief Core headers and types
 * 
 * \author    Youngik Yang
 * \version   0.2
 * \date      2013
 * \bug       Not known.
 * \warning   None.
 * \date      Modified on Fri 2013-12-13 03:27:04 PM
 * \copyright J. Craig Venter Institute.
 */

#ifndef __CORE_H__
#define __CORE_H__

/* core headers */
#include <iostream>
//#include <boost/lexical_cast.hpp>
//#include <boost/filesystem/path.hpp>
//#include <boost/filesystem/operations.hpp>
#include <tr1/unordered_map>
#include "param.h"
#include "sequence.h"
#include "gsa.h"
#include "path.h"
#include "alignsummary.h"
#include "filter.h"
#include "utils.h"
#include "semiglobal.h"
#include "default.h"
#include "meminfo.h"

enum DateType { AA_TYPE, DNA_TYPE };
enum OutputType { ALIGN, PROFILE, SEQUENCE, PLACE_TXT, PLACE_BIN };
const std::string SpaFiles[] = { "spa.align.aac", "spa.profile.aac", "spa.fasta", "spa.place.txt", "spa.place.bin" };

#endif
