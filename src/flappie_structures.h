/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#pragma once
#ifndef FLAPPIE_STRUCTURES_H
#define FLAPPIE_STRUCTURES_H

#include <stddef.h>
#include "flappie_matrix.h"

typedef struct {
    char * uuid;
    size_t n;
    size_t start;
    size_t end;
    float *raw;
} raw_table;

struct _raw_basecall_info {
    float score;
    raw_table rt;

    char *basecall;
    char *quality;
    size_t basecall_length;
    flappie_imatrix trace;

    int *pos;
    size_t nblock;
};

void free_raw_table(raw_table * tbl);
void free_raw_basecall_info(struct _raw_basecall_info * ptr);
#endif /* FLAPPIE_STRUCTURES_H */
