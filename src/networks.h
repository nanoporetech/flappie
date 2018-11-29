/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#pragma once
#ifndef NETWORKS_H
#    define NETWORKS_H
#    include <stdbool.h>
#    include "flappie_matrix.h"
#    include "flappie_structures.h"

typedef flappie_matrix (*transition_function_ptr)(const raw_table, float);

enum model_type {
    //FLAPPIE_MODEL_R94_PCR = 0,
    FLAPPIE_MODEL_R941_NATIVE = 0,
    FLAPPIE_MODEL_R941_5mC,
    FLAPPIE_MODEL_R10C_PCR,
    FLAPPIE_MODEL_INVALID};

static const enum model_type flappie_nmodel = FLAPPIE_MODEL_INVALID;
enum model_type get_flappie_model_type(const char *modelstr);
const char *flappie_model_string(const enum model_type model);
const char *flappie_model_description(const enum model_type model);
transition_function_ptr get_transition_function(const enum model_type model);

flappie_matrix flipflop_transitions(const raw_table signal, float temperature, enum model_type model);
flappie_matrix flipflop_transitions_r941native(const raw_table signal, float temperature);
flappie_matrix flipflop_transitions_r941native5mC(const raw_table signal, float temperature);
flappie_matrix flipflop_transitions_r10Cpcr(const raw_table signal, float temperature);


#endif    /* NETWORKS_H */
