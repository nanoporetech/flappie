/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#include <err.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "flappie_util.h"

int main(int argc, char * argv[]){
    if(4 != argc){
        fprintf(stderr, "Usage: read_crp filename nrow ncol\n");
        exit(EXIT_FAILURE);
    }
    flappie_matrix mat = read_flappie_matrix(argv[1]);
    if(NULL == mat){
        errx(EXIT_FAILURE, "Failed to read matrix from '%s'\n", argv[1]);
    }

    int nr = atoi(argv[2]);
    int nc = atoi(argv[3]);
    fprint_flappie_matrix(stdout, NULL, mat, nr, nc, false);

    return EXIT_SUCCESS;
}
