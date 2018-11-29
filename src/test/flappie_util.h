/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#ifndef FLAPPIE_MATRIX_UTIL
#    define FLAPPIE_MATRIX_UTIL

#    include <flappie_matrix.h>
#    include <stdio.h>

size_t write_flappie_matrix(const char * fn, const_flappie_matrix mat);
size_t write_flappie_matrix_to_handle(FILE * fh, const_flappie_matrix mat);
flappie_matrix read_flappie_matrix_from_handle(FILE * fh);
flappie_matrix read_flappie_matrix(char const * fn);

flappie_matrix random_flappie_matrix(size_t nr, size_t nc, float lower,
                                       float upper);

#endif                          /* FLAPPIE_MATRIX_UTIL */
