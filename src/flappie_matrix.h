/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#pragma once
#ifndef FLAPPIE_MATRIX_H
#    define FLAPPIE_MATRIX_H

#    include <immintrin.h>
#    include <stdbool.h>
#    include <stdint.h>
#    include <stdio.h>

typedef struct {
    size_t nr, nrq, nc, stride;
    union {
        __m128 *v;
        float *f;
    } data;
} _Mat;

typedef struct {
    size_t nr, nrq, nc, stride;
    union {
        __m128i *v;
        int32_t *f;
    } data;
} _iMat;

typedef _Mat *flappie_matrix;
typedef _iMat *flappie_imatrix;
typedef _Mat const *const_flappie_matrix;
typedef _iMat const *const_flappie_imatrix;

flappie_matrix make_flappie_matrix(size_t nr, size_t nc);
flappie_matrix remake_flappie_matrix(flappie_matrix M, size_t nr, size_t nc);
flappie_matrix copy_flappie_matrix(const_flappie_matrix mat);
flappie_matrix free_flappie_matrix(flappie_matrix mat);
void zero_flappie_matrix(flappie_matrix M);
flappie_matrix mat_from_array(const float *x, size_t nr, size_t nc);
float * array_from_flappie_matrix(const_flappie_matrix mat);
void fprint_flappie_matrix(FILE * fh, const char *header,
                            const_flappie_matrix mat, size_t nr, size_t nc,
                            bool include_padding);
bool equality_flappie_matrix(const_flappie_matrix mat1,
                              const_flappie_matrix mat2, const float tol);
bool validate_flappie_matrix(flappie_matrix mat, float lower,
                              const float upper, const float maskval,
                              const bool only_finite, const char *file,
                              const int line);

flappie_imatrix make_flappie_imatrix(size_t nr, size_t nc);
flappie_imatrix remake_flappie_imatrix(flappie_imatrix M, size_t nr, size_t nc);
flappie_imatrix copy_flappie_imatrix(const_flappie_imatrix mat);
flappie_imatrix free_flappie_imatrix(flappie_imatrix mat);
void zero_flappie_imatrix(flappie_imatrix M);

flappie_matrix affine_map(const_flappie_matrix X, const_flappie_matrix W,
                           const_flappie_matrix b, flappie_matrix C);
flappie_matrix affine_map2(const_flappie_matrix Xf, const_flappie_matrix Xb,
                            const_flappie_matrix Wf, const_flappie_matrix Wb,
                            const_flappie_matrix b, flappie_matrix C);
void row_normalise_inplace(flappie_matrix C);
void log_row_normalise_inplace(flappie_matrix C);

float min_flappie_matrix(const_flappie_matrix mat);
float max_flappie_matrix(const_flappie_matrix mat);

bool validate_ivector(int *vec, const size_t n, const int lower,
                      const int upper, const char *file, const int line);

bool validate_vector(float *vec, const size_t n, const float lower,
                     const float upper, const char *file, const int line);

void shift_scale_matrix_inplace(flappie_matrix sigmat, float shift, float scale);
void clip_matrix_inplace(flappie_matrix C, float thresh);
void filter_matrix_inplace(flappie_matrix C, float fill_val,float thresh);
void difference_matrix_inplace(flappie_matrix C, float val);


#endif                          /* FLAPPIE_MATRIX_H */
