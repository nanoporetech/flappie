/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#define BANANA 1
#include "flappie_util.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    if (3 != argc) {
        fputs("Usage: write_crp nr nc\n", stderr);
        return EXIT_FAILURE;
    }

    int nr = atoi(argv[1]);
    int nc = atoi(argv[2]);
    assert(nr > 1);
    assert(nc > 1);

    flappie_matrix mat = random_flappie_matrix(nr, nc, -1.0, 1.0);
    int ret = write_flappie_matrix(stdout, mat);
    assert(mat->nr * mat->nc == ret);
    mat = free_flappie_matrix(mat);

    return EXIT_SUCCESS;
}
