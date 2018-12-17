/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#ifdef __APPLE__
#    include <Accelerate/Accelerate.h>
#else
#    include <cblas.h>
#endif
#include <float.h>
#include <math.h>
#include "flappie_matrix.h"
#include "flappie_stdlib.h"
#include "util.h"

flappie_matrix make_flappie_matrix(size_t nr, size_t nc) {
    assert(nr > 0);
    assert(nc > 0);
    // Matrix padded so row length is multiple of 4
    const size_t nrq = (size_t)ceil(nr / 4.0);
    flappie_matrix mat = malloc(sizeof(*mat));
    RETURN_NULL_IF(NULL == mat, NULL);

    mat->nr = nr;
    mat->nrq = nrq;
    mat->nc = nc;
    mat->stride = nrq * 4;

    {
        // Check for overflow to please Coverity scanner
        size_t tmp1 = nrq * sizeof(__m128);
        size_t tmp2 = tmp1 * nc;
        if (tmp1 != 0 && tmp2 / tmp1 != nc) {
            // Have overflow in memory allocation
            free(mat);
            return NULL;
        }
    }

    if (0 != flappie_memalign((void **)&mat->data.v, 16, nrq * nc * sizeof(__m128))) {
        warnx("Error allocating memory in %s.\n", __func__);
        free(mat);
        return NULL;
    }
    memset(mat->data.v, 0, nrq * nc * sizeof(__m128));
    return mat;
}


flappie_matrix remake_flappie_matrix(flappie_matrix M, size_t nr, size_t nc) {
    // Could be made more efficient when there is sufficent memory already allocated
    if ((NULL == M) || (M->nr != nr) || (M->nc != nc)) {
        M = free_flappie_matrix(M);
        M = make_flappie_matrix(nr, nc);
    }
    return M;
}


flappie_matrix copy_flappie_matrix(const_flappie_matrix M){
    RETURN_NULL_IF(NULL == M, NULL);
    flappie_matrix C = make_flappie_matrix(M->nr, M->nc);
    RETURN_NULL_IF(NULL == C, NULL);
    memcpy(C->data.f, M->data.f, sizeof(__m128) * C->nrq * C->nc);
    return C;
}


void zero_flappie_matrix(flappie_matrix M) {
    if (NULL == M) {
        return;
    }
    memset(M->data.f, 0, M->stride * M->nc * sizeof(float));
}


flappie_matrix mat_from_array(const float *x, size_t nr, size_t nc) {
    flappie_matrix res = make_flappie_matrix(nr, nc);
    RETURN_NULL_IF(NULL == res, NULL);

    for (size_t col = 0; col < nc; col++) {
        memcpy(res->data.f + col * res->stride, x + col * nr,
               nr * sizeof(float));
    }
    return res;
}

float * array_from_flappie_matrix(const_flappie_matrix mat){
    RETURN_NULL_IF(NULL == mat, NULL);

    const size_t nelt = mat->nr * mat->nc;
    float * res = calloc(nelt, sizeof(float));
    RETURN_NULL_IF(NULL == res, NULL);

    for(size_t c=0 ; c < mat->nc ; c++){
        const size_t offset_out = c * mat->nr;
        const size_t offset_in = c * mat->stride;
        for(size_t r=0 ; r < mat->nr ; r++){
            res[offset_out + r] = mat->data.f[offset_in + r];
        }
    }

    return res;
}


void fprint_flappie_matrix(FILE * fh, const char *header,
                            const_flappie_matrix mat, size_t nr, size_t nc,
                            bool include_padding) {
    assert(NULL != fh);
    assert(NULL != mat);
    const size_t rlim = include_padding ? mat->stride : mat->nr;

    if (nr <= 0 || nr > rlim) {
        nr = rlim;
    }
    if (nc <= 0 || nc > mat->nc) {
        nc = mat->nc;
    }

    if (NULL != header) {
        int ret = fputs(header, fh);
        if (EOF == ret || ret < 0) {
            return;
        }
        fputc('\n', fh);
    }
    for (size_t c = 0; c < nc; c++) {
        const size_t offset = c * mat->stride;
        fprintf(fh, "%4zu : % 12e", c, mat->data.f[offset]);
        for (size_t r = 1; r < nr; r++) {
            fprintf(fh, "  % 12e", mat->data.f[offset + r]);
        }
        fputc('\n', fh);
    }
}

flappie_matrix free_flappie_matrix(flappie_matrix mat) {
    if (NULL != mat) {
        free(mat->data.v);
        free(mat);
    }
    return NULL;
}

bool validate_flappie_matrix(flappie_matrix mat, float lower,
                              const float upper, const float maskval,
                              const bool only_finite, const char *file,
                              const int line) {
#ifdef NDEBUG
    return true;
}
#else
    if (NULL == mat) {
        return false;
    }
    assert(NULL != mat->data.f);
    assert(mat->nc > 0);
    assert(mat->nr > 0);
    assert(mat->stride > 0 && mat->stride >= mat->nr);
    assert(mat->nrq * 4 == mat->stride);

    const size_t nc = mat->nc;
    const size_t nr = mat->nr;
    const size_t ld = mat->stride;

    //  Masked values correct
    if (!isnan(maskval)) {
        for (size_t c = 0; c < nc; ++c) {
            const size_t offset = c * ld;
            for (size_t r = nr; r < ld; ++r) {
                if (maskval != mat->data.f[offset + r]) {
                    warnx
                        ("%s:%d  Matrix entry [%zu,%zu] = %f violates masking rules\n",
                         file, line, r, c, mat->data.f[offset + r]);
                    return false;
                }
            }
        }
    }
    //  Check finite
    if (only_finite) {
        for (size_t c = 0; c < nc; ++c) {
            const size_t offset = c * ld;
            for (size_t r = 0; r < nr; ++r) {
                if (!isfinite(mat->data.f[offset + r])) {
                    warnx
                        ("%s:%d  Matrix entry [%zu,%zu] = %f contains a non-finite value\n",
                         file, line, r, c, mat->data.f[offset + r]);
                    return false;
                }
            }
        }
    }
    //  Lower bound
    if (!isnan(lower)) {
        for (size_t c = 0; c < nc; ++c) {
            const size_t offset = c * ld;
            for (size_t r = 0; r < nr; ++r) {
                if (mat->data.f[offset + r] + FLT_EPSILON < lower) {
                    warnx
                        ("%s:%d  Matrix entry [%zu,%zu] = %f (%e) violates lower bound\n",
                         file, line, r, c, mat->data.f[offset + r],
                         mat->data.f[offset + r] - lower);
                    return false;
                }
            }
        }
    }
    //  Upper bound
    if (!isnan(upper)) {
        for (size_t c = 0; c < nc; ++c) {
            const size_t offset = c * ld;
            for (size_t r = 0; r < nr; ++r) {
                if (mat->data.f[offset + r] > upper + FLT_EPSILON) {
                    warnx
                        ("%s:%d  Matrix entry [%zu,%zu] = %f (%e) violates upper bound\n",
                         file, line, r, c, mat->data.f[offset + r],
                         mat->data.f[offset + r] - upper);
                    return false;
                }
            }
        }
    }

    return true;
}
#endif /* NDEBUG */


/**  Check whether two matrices are equal within a given tolerance
 *
 *  @param mat1 A `flappie_matrix` to compare
 *  @param mat2 A `flappie_matrix` to compare
 *  @param tol Absolute tolerance to which elements of the matrix should agree
 *
 *  Notes:
 *    The tolerance is absolute; this may not be desirable in all circumstances.
 *    The convention used here is that of equality '=='.  The standard C
 *    sorting functions expect the convention of 0 being equal and non-equality
 *    being defined by negative (less than) and positive (greater than).
 *
 *  @return A boolean of whether the two matrices are equal.
 **/
bool equality_flappie_matrix(const_flappie_matrix mat1,
                              const_flappie_matrix mat2, const float tol) {
    if (NULL == mat1 || NULL == mat2) {
        // One or both matrices are NULL
        if (NULL == mat1 && NULL == mat2) {
            return true;
        }
        return false;
    }
    // Given non-NULL matrices, they should always contain data
    assert(NULL != mat1->data.f);
    assert(NULL != mat2->data.f);

    if (mat1->nc != mat2->nc || mat1->nr != mat2->nr) {
        // Dimension mismatch
        return false;
    }
    //  Given equal dimensions, the following should alway hold
    assert(mat1->nrq == mat2->nrq);

    for (size_t c = 0; c < mat1->nc; ++c) {
        const size_t offset = c * mat1->stride;
        for (size_t r = 0; r < mat1->nr; ++r) {
            if (fabsf(mat1->data.f[offset + r] - mat2->data.f[offset + r]) >
                tol) {
                return false;
            }
        }
    }

    return true;
}


flappie_imatrix make_flappie_imatrix(size_t nr, size_t nc) {
    assert(nr > 0);
    assert(nc > 0);
    // Matrix padded so row length is multiple of 4
    const size_t nrq = (size_t)ceil(nr / 4.0);
    flappie_imatrix mat = malloc(sizeof(*mat));
    RETURN_NULL_IF(NULL == mat, NULL);

    mat->nr = nr;
    mat->nrq = nrq;
    mat->nc = nc;
    mat->stride = nrq * 4;

    if (0 != flappie_memalign((void **)&mat->data.v, 16, nrq * nc * sizeof(__m128i))) {
        warnx("Error allocating memory in %s.\n", __func__);
        free(mat);
        return NULL;
    }
    memset(mat->data.v, 0, nrq * nc * sizeof(__m128));
    return mat;
}


flappie_imatrix remake_flappie_imatrix(flappie_imatrix M, size_t nr, size_t nc) {
    // Could be made more efficient when there is sufficent memory already allocated
    if ((NULL == M) || (M->nr != nr) || (M->nc != nc)) {
        M = free_flappie_imatrix(M);
        M = make_flappie_imatrix(nr, nc);
    }
    return M;
}


flappie_imatrix copy_flappie_imatrix(const_flappie_imatrix M){
    RETURN_NULL_IF(NULL == M, NULL);
    flappie_imatrix C = make_flappie_imatrix(M->nr, M->nc);
    RETURN_NULL_IF(NULL == C, NULL);
    memcpy(C->data.f, M->data.f, sizeof(__m128i) * C->nrq * C->nc);
    return C;
}


flappie_imatrix free_flappie_imatrix(flappie_imatrix mat) {
    if (NULL != mat) {
        free(mat->data.v);
        free(mat);
    }
    return NULL;
}


void zero_flappie_imatrix(flappie_imatrix M) {
    if (NULL == M) {
        return;
    }
    memset(M->data.f, 0, M->stride * M->nc * sizeof(int));
}


int32_t * array_from_flappie_imatrix(const_flappie_imatrix mat){
    RETURN_NULL_IF(NULL == mat, NULL);

    const size_t nelt = mat->nr * mat->nc;
    int32_t * res = calloc(nelt, sizeof(int32_t));
    RETURN_NULL_IF(NULL == res, NULL);

    for(size_t c=0 ; c < mat->nc ; c++){
        const size_t offset_out = c * mat->nr;
        const size_t offset_in = c * mat->stride;
        for(size_t r=0 ; r < mat->nr ; r++){
            res[offset_out + r] = mat->data.f[offset_in + r];
        }
    }

    return res;
}


flappie_matrix affine_map(const_flappie_matrix X, const_flappie_matrix W,
                           const_flappie_matrix b, flappie_matrix C) {
    /*  Affine transform C = W^t X + b
     *  X is [nr, nc]
     *  W is [nr, nk]
     *  b is [nk]
     *  C is [nk, nc] or NULL.  If NULL then C is allocated.
     */
    RETURN_NULL_IF(NULL == X, NULL);

    assert(NULL != W);
    assert(NULL != b);
    assert(W->nr == X->nr);

    C = remake_flappie_matrix(C, W->nc, X->nc);
    RETURN_NULL_IF(NULL == C, NULL);

    /* Copy bias */
    for (size_t c = 0; c < C->nc; c++) {
        memcpy(C->data.v + c * C->nrq, b->data.v, C->nrq * sizeof(__m128));
    }

    /* Affine transform */
    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, W->nc, X->nc, W->nr,
                1.0, W->data.f, W->stride, X->data.f, X->stride, 1.0,
                C->data.f, C->stride);

    return C;
}


flappie_matrix affine_map2(const_flappie_matrix Xf, const_flappie_matrix Xb,
                            const_flappie_matrix Wf, const_flappie_matrix Wb,
                            const_flappie_matrix b, flappie_matrix C) {
    RETURN_NULL_IF(NULL == Xf, NULL);
    RETURN_NULL_IF(NULL == Xb, NULL);

    assert(NULL != Wf);
    assert(NULL != Wb);
    assert(NULL != b);
    assert(Wf->nr == Xf->nr);
    assert(Wb->nr == Xb->nr);
    assert(Xf->nc == Xb->nc);
    assert(Wf->nc == Wb->nc);
    C = remake_flappie_matrix(C, Wf->nc, Xf->nc);
    RETURN_NULL_IF(NULL == C, NULL);

    /* Copy bias */
    for (size_t c = 0; c < C->nc; c++) {
        memcpy(C->data.v + c * C->nrq, b->data.v, C->nrq * sizeof(__m128));
    }

    /* Affine transform -- forwards */
    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, Wf->nc, Xf->nc, Wf->nr,
                1.0, Wf->data.f, Wf->stride, Xf->data.f, Xf->stride, 1.0,
                C->data.f, C->stride);
    /* Affine transform -- backwards */
    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, Wb->nc, Xb->nc, Wb->nr,
                1.0, Wb->data.f, Wb->stride, Xb->data.f, Xb->stride, 1.0,
                C->data.f, C->stride);
    return C;
}


void row_normalise_inplace(flappie_matrix C) {
    if (NULL == C) {
        // Input NULL due to earlier failure.  Propagate
        return;
    }
    const size_t i = C->stride - C->nr;
    const __m128 mask = _mm_cmpgt_ps(_mm_set_ps(i >= 1, i >= 2, i >= 3, 0), _mm_set1_ps(0.0f));
    for (size_t col = 0; col < C->nc; col++) {
        const size_t offset = col * C->nrq;
        __m128 sum = C->data.v[offset];
        for (size_t row = 1; row < C->nrq; row++) {
            sum += C->data.v[offset + row];
        }
        sum -= _mm_and_ps(C->data.v[offset + C->nrq - 1], mask);
        const __m128 psum = _mm_hadd_ps(sum, sum);
        const __m128 tsum = _mm_hadd_ps(psum, psum);

        const __m128 tsum_recip = _mm_set1_ps(1.0f) / tsum;
        for (size_t row = 0; row < C->nrq; row++) {
            C->data.v[offset + row] *= tsum_recip;
        }
    }
}


void log_row_normalise_inplace(flappie_matrix C){
    if(NULL == C){
        // Input NULL due to earlier failure.  Propagate
        return;
    }

    for (size_t col=0 ; col < C->nc; col++) {
        const size_t offset = col * C->stride;
        float row_logsum = C->data.f[offset];
        for(size_t row=1 ; row < C->nr ; row++){
            row_logsum = logsumexpf(row_logsum, C->data.f[offset + row]);
        }

        for(size_t row=0 ; row < C->nr ; row++){
            C->data.f[offset + row] -= row_logsum;
        }
    }
}


float max_flappie_matrix(const_flappie_matrix x) {
    if (NULL == x) {
        // Input NULL due to earlier failure.  Propagate
        return NAN;
    }
    float amax = x->data.f[0];
    for (size_t col = 0; col < x->nc; col++) {
        const size_t offset = col * x->stride;
        for (size_t r = 0; r < x->nr; r++) {
            if (amax < x->data.f[offset + r]) {
                amax = x->data.f[offset + r];
            }
        }
    }
    return amax;
}


float min_flappie_matrix(const_flappie_matrix x) {
    if (NULL == x) {
        // Input NULL due to earlier failure.  Propagate
        return NAN;
    }
    float amin = x->data.f[0];
    for (size_t col = 0; col < x->nc; col++) {
        const size_t offset = col * x->stride;
        for (size_t r = 0; r < x->nr; r++) {
            if (amin < x->data.f[offset + r]) {
                amin = x->data.f[offset + r];
            }
        }
    }
    return amin;
}


int argmax_flappie_matrix(const_flappie_matrix x) {
    if (NULL == x) {
        // Input NULL due to earlier failure.  Propagate
        return -1;
    }
    float amax = x->data.f[0];
    size_t imax = 0;

    for (size_t col = 0; col < x->nc; col++) {
        const size_t offset = col * x->stride;
        for (size_t r = 0; r < x->nr; r++) {
            if (amax < x->data.f[offset + r]) {
                amax = x->data.f[offset + r];
                imax = offset + r;
            }
        }
    }
    return imax;
}


int argmin_flappie_matrix(const_flappie_matrix x) {
    if (NULL == x) {
        // Input NULL due to earlier failure.  Propagate
        return -1;
    }
    float amin = x->data.f[0];
    size_t imin = 0;

    for (size_t col = 0; col < x->nc; col++) {
        const size_t offset = col * x->stride;
        for (size_t r = 0; r < x->nr; r++) {
            if (amin < x->data.f[offset + r]) {
                amin = x->data.f[offset + r];
                imin = offset + r;
            }
        }
    }
    return imin;
}


bool validate_vector(float *vec, const size_t n, const float lower,
                     const float upper, const char *file, const int line) {
#ifdef NDEBUG
    return true;
}
#else
    if (NULL == vec) {
        return false;
    }
    //  Lower bound
    if (!isnan(lower)) {
        for (size_t i = 0; i < n; ++i) {
            if (lower > vec[i]) {
                warnx("%s:%d  Vector entry %zu = %f violates lower bound\n",
                      file, line, i, vec[i]);
                return false;
            }
        }
    }
    //  Upper bound
    if (!isnan(upper)) {
        for (size_t i = 0; i < n; ++i) {
            if (upper < vec[i]) {
                warnx("%s:%d  Vector entry %zu = %f violates upper bound\n",
                      file, line, i, vec[i]);
                return false;
            }
        }
    }

    return true;
}
#endif /* NDEBUG */

bool validate_ivector(int *vec, const size_t n, const int lower, const int upper,
                      const char *file, const int line) {
#ifdef NDEBUG
    return true;
}
#else
    if (NULL == vec) {
        return false;
    }
    //  Lower bound
    for (size_t i = 0; i < n; ++i) {
        if (lower > vec[i]) {
            warnx("%s:%d  Vector entry %zu = %d violates lower bound\n", file,
                  line, i, vec[i]);
            return false;
        }
    }

    //  Upper bound
    for (size_t i = 0; i < n; ++i) {
        if (upper < vec[i]) {
            warnx("%s:%d  Vector entry %zu = %d violates upper bound\n", file,
                  line, i, vec[i]);
            return false;
        }
    }

    return true;
}
#endif /* NDEBUG */


/**  Shift-Scale a matrix elementwise
 *
 *   x[i] := (x[i] - shift) / scale
 *   Matrix updated in-place.
 *
 *   @param C     Matrix to transform [in/out]
 *   @param shift
 *   @param scale
 *
 *   @returns void
 **/
void shift_scale_matrix_inplace(flappie_matrix C, float shift, float scale){
    RETURN_NULL_IF(NULL == C,);
    for(size_t c=0 ; c < C->nc ; c++){
        const size_t offset = c * C->stride;
        for(size_t r=0 ; r < C->nr ; r++){
            C->data.f[offset + r] = (C->data.f[offset + r] - shift) / scale;
        }
    }
}


/**  Clip matrix into range
 *
 *   Clip elements of matrix into [-thresh, thresh].
 *   Matrix updated in-place.
 *
 *   @param x      Matrix to clip [in/out]
 *   @param thresh Threshold
 *
 *   @returns void
 **/
void clip_matrix_inplace(flappie_matrix C, float thresh){
    RETURN_NULL_IF(NULL == C,);
    for(size_t c=0 ; c < C->nc ; c++){
        const size_t offset = c * C->stride;
        for(size_t r=0 ; r < C->nr ; r++){
            const float obs = C->data.f[offset + r];
            const float val = fminf(thresh, fabsf(obs));
            C->data.f[offset + r] = copysign(val, obs);
        }
    }
}


/**  Filter absolutely large values from matrix
 *
 *   Replaces elements of matrix whose absolute value exceeds a threshhold.
 *   Matrix updated in-place
 *
 *   @param x        Matrix to filter [in/out]
 *   @param fill_val Value to replace filtered elements by
 *   @param thresh   Threshold
 *
 *   @returns void
 **/
void filter_matrix_inplace(flappie_matrix C, float fill_val, float thresh){
    RETURN_NULL_IF(NULL == C,);

    for(size_t c=0 ; c < C->nc ; c++){
        const size_t offset = c * C->stride;
        for(size_t r=0 ; r < C->nr ; r++){
            const float val = fabsf(C->data.f[offset + r]);
            if(val > thresh){
                C->data.f[offset + r] = fill_val;
            }
        }
    }
}


/**  Sliding differences across matrix
 *
 *   Calculated x[i] := x[i + 1] - x[i]
 *   Array updated in-place with final element set to `val`
 *
 *
 *   @param x Matrix to difference [in/out]
 *   @param val Value to pad with
 *
 *   @return  void
 **/
void difference_matrix_inplace(flappie_matrix C, float val){
    RETURN_NULL_IF(NULL == C,);

    for(size_t c=1 ; c < C->nc ; c++){
        const size_t offset = c * C->stride;
        const size_t offsetM1 = (c - 1) * C->stride;
        for(size_t r=0 ; r < C->nr ; r++){
            C->data.f[offsetM1 + r] = C->data.f[offset + r] - C->data.f[offsetM1 + r];
        }
    }
    {
        const size_t offset = (C->nc - 1) * C->stride;
        for(size_t r=0 ; r < C->nr ; r++){
            C->data.f[offset + r] = val;
        }
    }
}
