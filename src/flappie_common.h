#pragma once
#ifndef FLAPPIE_COMMON_H
#define FLAPPIE_COMMON_H

#include "flappie_structures.h"

raw_table trim_and_segment_raw(raw_table rt, size_t trim_start, size_t trim_end, size_t varseg_chunk, float varseg_thresh);
raw_table trim_raw_by_mad(raw_table rt, size_t chunk_size, float proportion);

#endif /* FLAPPIE_COMMON_H */
