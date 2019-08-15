/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#pragma once
#ifndef DECODE_H
#    define DECODE_H
#    include <stdbool.h>
#    include "flappie_matrix.h"
#    include "flappie_structures.h"

static const char base_lookup[5] = {'A', 'C', 'G', 'T', 'Z' };
static inline char basechar(int b){
    return base_lookup[b];
}

float argmax_decoder(const_flappie_matrix logpost, int *seq);
char * collapse_repeats(int const * path, size_t npos, int modbase);
size_t change_positions(int const * path, size_t npos, int * chpos);

float decode_crf_flipflop(const_flappie_matrix trans, bool combine_stays, int * path, float * qpath);
float decode_runlength(const_flappie_matrix param, int * path);
float decode_crf_runlength(const_flappie_matrix transparam, int * path);
float constrained_crf_flipflop(const_flappie_matrix post, int * path);

size_t runlengths_mean(const_flappie_matrix param, const int * path, int * runlength);
size_t runlengths_unit(const_flappie_matrix param, const int * path, int * runlength);
char * runlength_to_basecall(const int * path, const int * runlength, size_t nblk);

flappie_matrix posterior_crf_flipflop(const_flappie_matrix trans, bool return_log);
flappie_matrix transpost_crf_flipflop(const_flappie_matrix trans, bool return_log);
flappie_matrix posterior_runlength(const_flappie_matrix param);
flappie_matrix transpost_crf_runlength(const_flappie_matrix trans);
flappie_imatrix trace_from_posterior(flappie_matrix tpost);

#endif                          /* DECODE_H */
