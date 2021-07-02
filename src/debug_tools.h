/*
 * debug_tools.h
 *
 *  Created on: Sep 16, 2019
 *      Author: bprather
 */

#ifndef SRC_DEBUG_TOOLS_H_
#define SRC_DEBUG_TOOLS_H_

#include "decs.h"

void print_matrix(char *name, REAL g[NDIM][NDIM]);
void print_matrix_c(char *name, _Complex REAL g[NDIM][NDIM]);
void print_vector(char *name, REAL v[NDIM]);


void check_ortho(REAL Econ[NDIM][NDIM], REAL Ecov[NDIM][NDIM]);
void check_u(REAL Ucon[NDIM], REAL Ucov[NDIM]);

/*
 Check that the coherency tensor N satisfies certain
 basic properties:
 k . N = N . k = 0
 hermitian
 evaluate the invariants: I, Q^2 + U^2, V^2
 */
void check_N(_Complex REAL N[NDIM][NDIM], REAL Kcon[NDIM],
    REAL gcov[NDIM][NDIM]);

#endif /* SRC_DEBUG_TOOLS_H_ */
