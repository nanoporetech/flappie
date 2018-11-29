/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#include "flappie_common.h"
#include "flappie_stdlib.h"
#include "util.h"

raw_table trim_and_segment_raw(raw_table rt, size_t trim_start, size_t trim_end, size_t varseg_chunk, float varseg_thresh) {
    RETURN_NULL_IF(NULL == rt.raw, (raw_table){0});

    rt = trim_raw_by_mad(rt, varseg_chunk, varseg_thresh);
    RETURN_NULL_IF(NULL == rt.raw, (raw_table){0});

    rt.start = (rt.n - rt.start) > trim_start ? rt.start + trim_start : rt.n;
    rt.end = (rt.end > trim_end) ? rt.end - trim_end : 0;

    if (rt.start >= rt.end) {
        free(rt.raw);
        return (raw_table){0};
    }

    return rt;
}

/**  Simple segmentation of a raw read by thresholding the MAD
 *
 *  The MAD of the raw signal is calculated for non-overlapping chunks and then
 *  thresholded to find regions at the beginning and end of the signal that have
 *  unusually low variation (generally a stall or open pore).  The threshhold is
 *  derived from the distribution of the calaculated MADs.
 *
 *  The threshold is chosen to be high since a single chunk above it will trigger
 *  the end of the trimming: the threshhold is chosen so it is unlikely to be
 *  exceeded in the leader but commonly exceeded in the main read.
 *
 *  @param rt Structure containing raw signal
 *  @param chunk_size Size of non-overlapping chunks
 *  @param perc  The quantile to be calculated to use for threshholding
 *
 *  @return A range structure containing new start and end for read
 **/
raw_table trim_raw_by_mad(raw_table rt, size_t chunk_size, float perc) {
    assert(chunk_size > 1);
    assert(perc >= 0.0 && perc <= 1.0);

    const size_t nsample = rt.end - rt.start;
    const size_t nchunk = nsample / chunk_size;
    // Truncation of end to be consistent with Sloika
    rt.end = nchunk * chunk_size;

    float *madarr = malloc(nchunk * sizeof(float));
    RETURN_NULL_IF(NULL == madarr, (raw_table){0});
    for (size_t i = 0; i < nchunk; i++) {
        madarr[i] = madf(rt.raw + rt.start + i * chunk_size, chunk_size, NULL);
    }
    quantilef(madarr, nchunk, &perc, 1);

    const float thresh = perc;
    for (size_t i = 0; i < nchunk; i++) {
        if (madarr[i] > thresh) {
            break;
        }
        rt.start += chunk_size;
    }
    for (size_t i = nchunk; i > 0; i--) {
        if (madarr[i - 1] > thresh) {
            break;
        }
        rt.end -= chunk_size;
    }
    assert(rt.end > rt.start);

    free(madarr);

    return rt;
}
