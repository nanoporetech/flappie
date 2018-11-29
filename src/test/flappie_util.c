/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#define BANANA 1
#define _BSD_SOURCE

#include <math.h>

#include "flappie_stdlib.h"
#include "flappie_util.h"

#define CRP_MAX_ROW (1<<16)
#define CRP_MAX_COL (1<<16)

/**  Simple writer for flappie_matrix structures
 *
 *  Writes a Scrappie matrix to a file point using hexadecimal ecoded
 *  floats.
 *
 *  @param fh A FILE pointer to write to.
 *  @param mat A matrix to write out.
 *
 *  @returns Number of elements written
 **/
size_t write_flappie_matrix_to_handle(FILE * fh, const_flappie_matrix mat) {
    if (NULL == fh || NULL == mat) {
        return 0;
    }

    fprintf(fh, "%zu\t%zu\n", mat->nr, mat->nc);

    size_t nelt = 0;
    for (size_t c = 0; c < mat->nc; ++c) {
        const size_t offset = c * mat->stride;
        int ret = fprintf(fh, "%a", mat->data.f[offset + 0]);
        if (ret > 0) {
            nelt += 1;
        }
        for (size_t r = 1; r < mat->nr; ++r) {
            ret = fprintf(fh, "\t%a", mat->data.f[offset + r]);
            if (ret > 0) {
                nelt += 1;
            }
        }
        fputc('\n', fh);
    }

    assert(nelt == mat->nr * mat->nc);
    return nelt;
}

size_t write_flappie_matrix(const char * fn, const_flappie_matrix mat) {
    if(NULL == fn){
        return 0;
    }

    FILE * fh = fopen(fn, "w");
    if(NULL == fh){
        return 0;
    }

    size_t res = write_flappie_matrix_to_handle(fh, mat);
    fclose(fh);

    return res;
}

/**   Simple reader for flappie_matrix structures
 *
 *  @param fh A FILE pointer to read from.
 *
 *  @returns Matrix read or NULL on failure
 **/
flappie_matrix read_flappie_matrix_from_handle(FILE * fh) {
    RETURN_NULL_IF(NULL == fh, NULL);

    int inr, inc;
    int ret = fscanf(fh, "%d\t%d\n", &inr, &inc);
    if (2 != ret) {
        warnx("Invalid header line.\n");
        return NULL;
    }
    if (inr <= 0 || inc <= 0) {
        warnx("Number of rows or columns invalid (got %d %d).\n", inr, inc);
        return NULL;
    }
    if (inr > CRP_MAX_ROW || inc > CRP_MAX_COL) {
        warnx("Number of rows or columns too large (got %d %d).\n", inr, inc);
        return NULL;
    }

    size_t nr = inr;
    size_t nc = inc;

    flappie_matrix mat = make_flappie_matrix(nr, nc);
    if (NULL == mat) {
        warnx("Failed to allocate enough memory for matrix.\n");
        return NULL;
    }

    size_t nelt = 0;
    for (size_t c = 0; c < nc; ++c) {
        const size_t offset = c * mat->stride;
        ret = fscanf(fh, "%a", &mat->data.f[offset]);
        if (ret > 0) {
            nelt += 1;
        }
        for (size_t r = 1; r < nr; ++r) {
            ret = fscanf(fh, "\t%a", &mat->data.f[offset + r]);
            if (ret > 0) {
                nelt += 1;
            }
        }
    }

    if (nelt != nr * nc) {
        warnx
            ("Read incorrect number of elements. Got %zu but expecteding %zu x %zu\n",
             nelt, nr, nc);
        mat = free_flappie_matrix(mat);
        return NULL;
    }

    assert(validate_flappie_matrix
           (mat, NAN, NAN, 0.0, true, __FILE__, __LINE__));
    return mat;
}

flappie_matrix read_flappie_matrix(char const * fn) {
    RETURN_NULL_IF(NULL == fn, NULL);
    FILE *fh = fopen(fn, "r");
    RETURN_NULL_IF(NULL == fh, NULL);

    flappie_matrix mat = read_flappie_matrix_from_handle(fh);
    fclose(fh);
    return mat;
}

/**   Uniform number distributed [lower, upper]
 *
 *  @param lower Lower bound
 *  @param upper Upper bound
 *
 *  @returns A unformly distributed random number
 **/
float flappie_random_uniform(float lower, float upper) {
    // random() is awful but we don't require good random numbers
    return lower + (upper - lower) * (float)((double)random() / RAND_MAX);
}

/**  Matrix filed with random values
 *
 *   Creates new matrix whose elemenrs are picked from a uniform distribution
 *   on [lower, upper].
 *
 *  @param nr Number of rows of matrix
 *  @param nc Number of columns of matrix
 *  @param lower Lower bound for random values
 *  @param upper Upper bound for random values
 *
 *  @returns Matrix filled with random data, or NULL on error
 **/
flappie_matrix random_flappie_matrix(size_t nr, size_t nc, float lower, float upper) {
    assert(nr > 0);
    assert(nc > 0);
    assert(upper >= lower);

    flappie_matrix mat = make_flappie_matrix(nr, nc);
    RETURN_NULL_IF(NULL == mat, NULL);

    for (size_t c = 0; c < mat->nc; ++c) {
        const size_t offset = c * mat->stride;
        for (size_t r = 0; r < mat->nr; ++r) {
            mat->data.f[offset + r] = flappie_random_uniform(lower, upper);
        }
    }

    assert(validate_flappie_matrix
           (mat, lower, upper, 0.0, true, __FILE__, __LINE__));
    return mat;
}
