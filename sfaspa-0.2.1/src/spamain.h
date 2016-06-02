/** 
 * \file      spamain.h
 * \brief     Main program of SPA assembler
 * \author    Youngik Yang
 * \version   0.2
 * \date      2013
 * \date      Fri 2013-12-13 06:36:29 PM
 * \pre       Run prespa to generate graph input and suffix arrays.
 * \bug       None.
 * \warning   It may need large RAMs depending on data set.
 * \copyright J. Craig Venter Institute.
 */

#ifndef __SPA_MAIN_H__
#define __SPA_MAIN_H__

#include "cmdargs_simple.h"
#include "assembler.h"

/** 
 * Read assembler parameters
 */
void readParams(int argc, char **argv, Param &param);

#endif

