#pragma once
#ifndef FLAPPIE_STRUCTURES_H
#define FLAPPIE_STRUCTURES_H

#include <stddef.h>

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

    int *pos;
    size_t nblock;
};


#endif /* FLAPPIE_STRUCTURES_H */
