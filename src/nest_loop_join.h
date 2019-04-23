/**
 * @file    nest_loop_join.h
 * @author  Cagri Balkesen <cagri.balkesen@inf.ethz.ch>
 * @date    Sun Feb  5 20:12:56 2012
 * @version $Id: nest_loop_join.h 3017 2019-4-07 10:56:20Z bcagri $
 * 
 * @brief  The interface of nest loop join algorithm.
 *
 * (c) 2019, Liu Chaoyang, Chai Lab 
 *
 */

#ifndef _NEST_LOOP_JOIN_H
#define _NEST_LOOP_JOIN_H

#include "types.h" /* relation_t */

int64_t 
NLJ(relation_t *relR, relation_t *relS, int nthreads);

int64_t 
SNLJ(relation_t *relR, relation_t *relS, int nthreads);

int64_t 
SPNLJ(relation_t *relR, relation_t *relS, int nthreads);
#endif
