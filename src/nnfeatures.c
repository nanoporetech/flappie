/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#include <math.h>
#include <stdio.h>
#include "nnfeatures.h"
#include "flappie_stdlib.h"
#include "util.h"

flappie_matrix features_from_raw(const raw_table signal) {
    RETURN_NULL_IF(0 == signal.n, NULL);
    RETURN_NULL_IF(NULL == signal.raw, NULL);
    const size_t nsample = signal.end - signal.start;
    flappie_matrix sigmat = make_flappie_matrix(1, nsample);
    RETURN_NULL_IF(NULL == sigmat, NULL);

    const size_t offset = signal.start;
    for (size_t i = 0 ; i < nsample ; i++) {
        // Copy with stride 4 because of required padding for matrix
        sigmat->data.f[i * 4] = signal.raw[i + offset];
    }
    return sigmat;
}
