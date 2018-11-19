#pragma once
#ifndef DECODE_H
#    define DECODE_H
#    include <stdbool.h>
#    include "flappie_matrix.h"
#    include "flappie_structures.h"


float argmax_decoder(const_flappie_matrix logpost, int *seq);
char * collapse_repeats(int const * path, size_t npos, int modbase);

float decode_crf_flipflop(const_flappie_matrix trans, bool combine_stays, int * path);
float constrained_crf_flipflop(const_flappie_matrix post, int * path);

flappie_matrix posterior_crf_flipflop(const_flappie_matrix trans, bool return_log);
flappie_matrix transpost_crf_flipflop(const_flappie_matrix trans, bool return_log);

#endif                          /* DECODE_H */
