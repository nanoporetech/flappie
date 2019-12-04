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
#include <math.h>
#include "layers.h"
#include "flappie_stdlib.h"
#include "util.h"


/**  Apply tanh to a matrix element-wise
 *  @param C Matrix
 *
 **/
void tanh_activation_inplace(flappie_matrix C) {
    RETURN_NULL_IF(NULL == C, );
    for (size_t c = 0; c < C->nc; ++c) {
        const size_t offset = c * C->nrq;
        for (size_t r = 0; r < C->nrq; ++r) {
            C->data.v[offset + r] = TANHFV(C->data.v[offset + r]);
        }
    }
    (void)validate_flappie_matrix(C, -1.0, 1.0, 0.0, true, __FILE__, __LINE__);
}


/**  Apply exp to a matrix element-wise
 *  @param C Matrix
 *
 **/
void exp_activation_inplace(flappie_matrix C) {
    RETURN_NULL_IF(NULL == C, );
    for (size_t c = 0; c < C->nc; ++c) {
        const size_t offset = c * C->nrq;
        for (size_t r = 0; r < C->nrq; ++r) {
            C->data.v[offset + r] = EXPFV(C->data.v[offset + r]);
        }
    }
    (void)validate_flappie_matrix(C, 0.0, INFINITY, 1.0, true, __FILE__,
                                   __LINE__);
}


/**  Apply log to a matrix element-wise
 *  @param C Matrix
 *
 **/
void log_activation_inplace(flappie_matrix C) {
    RETURN_NULL_IF(NULL == C, );
    for (size_t c = 0; c < C->nc; ++c) {
        const size_t offset = c * C->nrq;
        for (size_t r = 0; r < C->nrq; ++r) {
            C->data.v[offset + r] = LOGFV(C->data.v[offset + r]);
        }
    }
}


/**  Apply ELU activation function to a matrix element-wise
 *  @param C Matrix
 *
 **/
void elu_activation_inplace(flappie_matrix C) {
    RETURN_NULL_IF(NULL == C, );
    for (size_t c = 0; c < C->nc; ++c) {
        const size_t offset = c * C->nrq;
        for (size_t r = 0; r < C->nrq; ++r) {
            C->data.v[offset + r] = ELUFV(C->data.v[offset + r]);
        }
    }
}


/** Apply robost log activation
 *
 *  Applies log(min_prob / nrow + (1 - min_prob) * x) elementwise to matrix
 *  where x in element and nrow is the number of rows
 *
 *  @param C Matrix
 *  @param min_prob  Minimum probability
 *
 **/
void robustlog_activation_inplace(flappie_matrix C, float min_prob) {
    assert(min_prob >= 0.0);
    assert(min_prob <= 1.0);
    RETURN_NULL_IF(NULL == C, );

    const size_t nblock = C->nc;
    const __m128 mpv = _mm_set1_ps(min_prob);
    const __m128 mpvm1 = _mm_set1_ps(1.0f - min_prob);
    for (size_t i = 0; i < nblock; i++) {
        const size_t offset = i * C->nrq;
        for (size_t r = 0; r < C->nrq; r++) {
            C->data.v[offset + r] =
                LOGFV(mpv + mpvm1 * C->data.v[offset + r]);
        }
    }
}


flappie_matrix embedding(int const * index, size_t n, const_flappie_matrix E, flappie_matrix C){
    RETURN_NULL_IF(NULL == index, NULL);

    const size_t nr = E->nr;
    const size_t nrq = E->nrq;
    const size_t nc = n;
    C = remake_flappie_matrix(C, nr, nc);
    RETURN_NULL_IF(NULL == C, NULL);

    for(size_t c=0 ; c < nc ; c++){
        assert(index[c] >= 0 && index[c] < E->nc);
        const size_t offsetC = c * nrq;
        const size_t offsetE = index[c] * nrq;
        for(size_t r=0 ; r < nrq ; r++){
            C->data.v[offsetC + r] = E->data.v[offsetE + r];
        }
    }

    return C;
}


flappie_matrix window(const_flappie_matrix input, size_t w, size_t stride) {
    RETURN_NULL_IF(NULL == input, NULL);
    assert(w > 0);
    const size_t wh = (w + 1) / 2;

    flappie_matrix output = make_flappie_matrix(input->nr * w,
                                                  (size_t)ceilf(input->nc /
                                                             (float)stride));
    RETURN_NULL_IF(NULL == output, NULL);

    for (size_t col = 0; col < output->nc; col++) {
        // First and last columns are special cases
        const size_t out_offset = col * output->stride;
        const int icol = (int)(col * stride);
        for (int i = 0, w1 = (icol - wh + 1); w1 <= icol + wh; w1++) {
            if (w1 < 0 || w1 >= input->nc) {
                i += input->nr;
                continue;
            }
            const size_t in_offset = w1 * input->stride;
            for (size_t row = 0; row < input->nr; row++, i++) {
                output->data.f[out_offset + i] = input->data.f[in_offset + row];
            }
        }
    }

    return output;
}


/**  Convolution of the input data
 *  @param X Input data matrix (features x nobs)
 *  @param W Filter matrix (winlen * features x nfilter)
 *
 *  The input is padded with zeros such that the resultant matrix has the
 *  same size as the input (under a stride of 1).
 *
 *  Note: The rows of the input matrix X are padded with zeros to make them
 *  a multiple of the SSE vector size (4).  The filter matrix must have been
 *  expanded accordingly.
 **/
flappie_matrix convolution(const_flappie_matrix X, const_flappie_matrix W,
                            const_flappie_matrix b, size_t stride,
                            flappie_matrix C) {
    RETURN_NULL_IF(NULL == X, NULL);
    assert(NULL != W);
    assert(NULL != b);
    assert(W->nc == b->nr);
    assert(stride > 0);
    // Window length of filter
    assert((W->nrq % X->nrq) == 0);
    const size_t winlen = W->nrq / X->nrq;
    const size_t nfilter = W->nc;
    // Padding -- right-hand side is longer when asymmetric padding is required
    const size_t padL = (winlen - 1) / 2;
    const size_t padR = winlen / 2;
    const size_t ncolC = iceil(X->nc, stride);
    C = remake_flappie_matrix(C, nfilter, ncolC);
    RETURN_NULL_IF(NULL == C, NULL);

    // Matrix strides
    const size_t ldC = C->stride;
    const size_t ldW = W->stride;
    const size_t ldX = X->stride;
    const size_t ldFeature = ldX;

    // Copy bias into result matrix
    for (size_t i = 0; i < C->nc; i++) {
        memcpy(C->data.v + i * C->nrq, b->data.v, C->nrq * sizeof(__m128));
    }

    // Left-hand side edge case where only part of the filter covers the input
    for (size_t w = 0; w < padL; w += stride) {
        const size_t offsetW = ldFeature * (padL - w);
        const size_t ncol = w / stride;
        cblas_sgemv(CblasColMajor, CblasTrans, W->nr - offsetW, W->nc,
                    1.0, W->data.f + offsetW, ldW,
                    X->data.f, 1, 1.0, C->data.f + ldC * ncol, 1);
    }

    // Number of columns of X already filled * ldC
    const size_t ncolsL_complete = iceil(padL, stride);
    const size_t offsetC_L = ldC * ncolsL_complete;
    // Because of stride, first valid filter may not start at beginning of X
    //const int shiftX_L = stride - (padL % stride);
    const size_t shiftX_L = ncolsL_complete * stride - padL;
    const size_t offsetX_L = shiftX_L * ldX;
    // Find multiple of stride greater or equal to winlen
    const size_t nstepC = iceil(winlen, stride);
    const size_t nstepX = stride * nstepC;

    for (size_t w = 0; w < winlen; w += stride) {
        //  Multiply reshaped X matrix by filter matrix
        //  The rows of X are padded by zeros to make a multiple of 4.
        //  Input matrix 'X'
        //   - stride is ldX * nstepX
        //   - offset by ldX * w (w cols)
        //   - Ncolumns is (X->nc - w) / nstepX + adjustment if a final window fits
        //  Filter matrix needs to be padded appropriately for the padding of X.
        //
        const size_t ncol_processed = ifloor(X->nc - shiftX_L - w, nstepX);
        const size_t initial_col = ifloor(w, stride);
        cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, W->nc,
                    ncol_processed, W->nr, 1.0, W->data.f, ldW,
                    X->data.f + ldX * w + offsetX_L, ldX * nstepX, 1.0,
                    C->data.f + ldC * initial_col + offsetC_L, ldC * nstepC);
    }

    // Right-hand side edge case where only part of the filter covers the input
    const size_t maxCol_reshape = ifloor(X->nc - shiftX_L, nstepX);
    const size_t remainder_reshape = (X->nc - shiftX_L) % nstepX;
    const size_t offsetC_R =
        offsetC_L + ldC * nstepC * (maxCol_reshape - 1) +
        ldC * (remainder_reshape / stride) + ldC;
    const size_t offsetX_R = (X->nc - winlen + 1) * ldX;
    // How far into padding is first block
    const int startR = stride - (padL + X->nc - winlen) % stride - 1;
    for (size_t w = startR; w < padR; w += stride) {
        const size_t offsetW = ldFeature * (w + 1);
        cblas_sgemv(CblasColMajor, CblasTrans, W->nr - offsetW, W->nc, 1.0,
                    W->data.f, ldW,
                    X->data.f + offsetX_R + ldX * w, 1, 1.0,
                    C->data.f + offsetC_R + ldC * (w / stride), 1);
    }

    assert(validate_flappie_matrix
           (C, NAN, NAN, 0.0, true, __FILE__, __LINE__));
    return C;
}


flappie_matrix feedforward_linear(const_flappie_matrix X,
                                   const_flappie_matrix W,
                                   const_flappie_matrix b, flappie_matrix C) {
    return affine_map(X, W, b, C);
}


flappie_matrix feedforward_tanh(const_flappie_matrix X,
                                 const_flappie_matrix W,
                                 const_flappie_matrix b, flappie_matrix C) {
    C = affine_map(X, W, b, C);
    RETURN_NULL_IF(NULL == C, NULL);

    tanh_activation_inplace(C);

    assert(validate_flappie_matrix
           (C, -1.0, 1.0, 0.0, true, __FILE__, __LINE__));
    return C;
}


flappie_matrix feedforward_exp(const_flappie_matrix X,
                                const_flappie_matrix W,
                                const_flappie_matrix b, flappie_matrix C) {
    C = affine_map(X, W, b, C);
    RETURN_NULL_IF(NULL == C, NULL);

    exp_activation_inplace(C);

    assert(validate_flappie_matrix
           (C, 0.0, NAN, 1.0, true, __FILE__, __LINE__));
    return C;
}


flappie_matrix residual(const_flappie_matrix X, const_flappie_matrix fX, flappie_matrix C) {
    RETURN_NULL_IF(NULL == X, NULL);
    RETURN_NULL_IF(NULL == fX, NULL);
    const size_t nr = X->nr;
    const size_t nrq = X->nrq;
    const size_t nc = X->nc;
    assert(nr == fX->nr);
    assert(nrq == fX->nrq);
    assert(nc == fX->nc);

    C = remake_flappie_matrix(C, nr, nc);
    RETURN_NULL_IF(NULL == C, NULL);

    for(size_t c=0 ; c < nc ; c++){
        const size_t offset = c * nrq;
        for(size_t r=0 ; r < nrq ; r++){
            C->data.v[offset + r] = X->data.v[offset + r] + fX->data.v[offset + r];
        }
    }

    return C;
}


void residual_inplace(const_flappie_matrix X, flappie_matrix fX) {
    RETURN_NULL_IF(NULL == X, );
    RETURN_NULL_IF(NULL == fX, );

    const size_t nrq = X->nrq;
    const size_t nc = X->nc;
    assert(X->nr == fX->nr);
    assert(nrq == fX->nrq);
    assert(nc == fX->nc);

    for(size_t c=0 ; c < nc ; c++){
        const size_t offset = c * nrq;
        for(size_t r=0 ; r < nrq ; r++){
            fX->data.v[offset + r] += X->data.v[offset + r];
        }
    }
}


flappie_matrix softmax(const_flappie_matrix X, const_flappie_matrix W,
                        const_flappie_matrix b, flappie_matrix C) {
    C = feedforward_exp(X, W, b, C);
    RETURN_NULL_IF(NULL == C, NULL);

    row_normalise_inplace(C);

    assert(validate_flappie_matrix
           (C, 0.0, 1.0, NAN, true, __FILE__, __LINE__));
    return C;
}


/**   Softmax with separate temperatures on weights and bias
 *
 *    Calculates softmax( A x / tempW + b / tempb ) as
 *    softmax( (A (x * tempb / tempW ) + b) / tempb )
 *
 *    @returns matrix containing softmax
 **/
flappie_matrix softmax_with_temperature(flappie_matrix X, const_flappie_matrix W,
                                         const_flappie_matrix b, float tempW, float tempb,
                                         flappie_matrix C) {
    RETURN_NULL_IF(NULL == X, NULL);

    shift_scale_matrix_inplace(X, 0.0f, tempW / tempb);

    C = feedforward_linear(X, W, b, C);
    RETURN_NULL_IF(NULL == C, NULL);

    shift_scale_matrix_inplace(C, 0.0f, tempb);
    exp_activation_inplace(C);
    row_normalise_inplace(C);

    assert(validate_flappie_matrix
           (C, 0.0, 1.0, NAN, true, __FILE__, __LINE__));
    return C;
}


flappie_matrix feedforward2_tanh(const_flappie_matrix Xf,
                                  const_flappie_matrix Xb,
                                  const_flappie_matrix Wf,
                                  const_flappie_matrix Wb,
                                  const_flappie_matrix b, flappie_matrix C) {
    C = affine_map2(Xf, Xb, Wf, Wb, b, C);
    RETURN_NULL_IF(NULL == C, NULL);

    tanh_activation_inplace(C);

    assert(validate_flappie_matrix(C, -1.0, 1.0, 0.0, true, __FILE__, __LINE__));
    return C;
}


flappie_matrix gru_forward(const_flappie_matrix X, const_flappie_matrix sW,
                            const_flappie_matrix sW2, flappie_matrix ostate) {
    RETURN_NULL_IF(NULL == X, NULL);

    assert(NULL != sW);
    assert(NULL != sW2);

    const size_t bsize = X->nc;
    const size_t size = sW2->nc;
    assert(X->nr == 3 * size);
    assert(sW->nr == size);
    assert(sW2->nr == size);
    assert(sW->nc == 2 * size);
    assert(sW2->nc == size);

    ostate = remake_flappie_matrix(ostate, size, bsize);
    RETURN_NULL_IF(NULL == ostate, NULL);

    flappie_matrix tmp = make_flappie_matrix(3 * size, 1);
    if(NULL == tmp){
        //  Memory allocation falled, clean-up and return
        free(ostate);
        return NULL;
    }

    /* First step state is zero.  Set second column of ostate to zero and use that */
    _Mat xCol, sCol1, sCol2;
    memset(ostate->data.v + ostate->nrq, 0, ostate->nrq * sizeof(__m128));
    xCol = *X;
    sCol1 = *ostate;
    sCol2 = *ostate;
    xCol.nc = sCol1.nc = sCol2.nc = 1;
    sCol1.data.v = ostate->data.v + ostate->nrq;
    sCol2.data.v = ostate->data.v;
    gru_step(&xCol, &sCol1, sW, sW2, tmp, &sCol2);
    for (size_t i = 1; i < bsize; i++) {
        xCol.data.v = X->data.v + i * X->nrq;
        sCol1.data.v = ostate->data.v + (i - 1) * ostate->nrq;
        sCol2.data.v = ostate->data.v + i * ostate->nrq;
        gru_step(&xCol, &sCol1, sW, sW2, tmp, &sCol2);
    }

    tmp = free_flappie_matrix(tmp);

    assert(validate_flappie_matrix
           (ostate, -1.0, 1.0, 0.0, true, __FILE__, __LINE__));
    return ostate;
}


flappie_matrix gru_backward(const_flappie_matrix X, const_flappie_matrix sW,
                             const_flappie_matrix sW2, flappie_matrix ostate) {
    RETURN_NULL_IF(NULL == X, NULL);
    assert(NULL != sW);
    assert(NULL != sW2);

    const size_t size = sW2->nc;
    const size_t bsize = X->nc;
    assert(X->nr == 3 * size);
    assert(sW->nr == size);
    assert(sW2->nr == size);
    assert(sW->nc == 2 * size);
    assert(sW2->nc == size);

    ostate = remake_flappie_matrix(ostate, size, bsize);
    RETURN_NULL_IF(NULL == ostate, NULL);

    flappie_matrix tmp = make_flappie_matrix(3 * size, 1);
    if(NULL == tmp){
        //  Memory allocation falled, clean-up and return
        free(ostate);
        return NULL;
    }

    /* First step state is zero.  Set first column of ostate to zero and use that */
    _Mat xCol, sCol1, sCol2;
    memset(ostate->data.v, 0, ostate->nrq * sizeof(__m128));
    xCol = *X;
    sCol1 = *ostate;
    sCol2 = *ostate;
    xCol.nc = sCol1.nc = sCol2.nc = 1;
    xCol.data.v = X->data.v + (X->nc - 1) * X->nrq;
    sCol1.data.v = ostate->data.v;
    sCol2.data.v = ostate->data.v + (ostate->nc - 1) * ostate->nrq;
    gru_step(&xCol, &sCol1, sW, sW2, tmp, &sCol2);
    for (size_t i = 1; i < bsize; i++) {
        const size_t index = bsize - i - 1;
        xCol.data.v = X->data.v + index * X->nrq;
        sCol1.data.v = ostate->data.v + (index + 1) * ostate->nrq;
        sCol2.data.v = ostate->data.v + index * ostate->nrq;
        gru_step(&xCol, &sCol1, sW, sW2, tmp, &sCol2);
    }

    tmp = free_flappie_matrix(tmp);

    assert(validate_flappie_matrix
           (ostate, -1.0, 1.0, 0.0, true, __FILE__, __LINE__));
    return ostate;
}


void gru_step(const_flappie_matrix x, const_flappie_matrix istate,
              const_flappie_matrix sW, const_flappie_matrix sW2,
              flappie_matrix xF, flappie_matrix ostate) {
    /* Perform a single GRU step
     * x      is [isize]
     * istate is [size]
     * xW     is [isize, 3 * size]
     * sW     is [size, 2 * size]
     * sW2    is [size, size]
     * bias   is [3 * size]
     * xF     is [3 * size]
     * ostate is [size]
     */
    assert(NULL != x);
    assert(NULL != sW);
    assert(NULL != sW2);
    const size_t size = istate->nr;
    assert(x->nr == 3 * size);
    assert(size % 4 == 0);  // Vectorisation assumes size divisible by 4
    const size_t sizeq = size / 4;
    assert(size == sW->nr);
    assert(2 * size == sW->nc);
    assert(size == sW2->nr);
    assert(size == sW2->nc);
    assert(3 * size == xF->nr);
    assert(size == ostate->nr);


    // Copy input vector = iW x + b to temporary vector
    memcpy(xF->data.v, x->data.v, x->nrq * sizeof(__m128));
    /*  Add sW * istate to first 2 * size elts of xF
     *  then apply gate function to get r and z
     */
    cblas_sgemv(CblasColMajor, CblasTrans, sW->nr, sW->nc, 1.0, sW->data.f,
                sW->stride, istate->data.f, 1, 1.0, xF->data.f, 1);
    for (size_t i = 0; i < (sizeq +sizeq); i++) {
        xF->data.v[i] = LOGISTICFV(xF->data.v[i]);
    }

    const __m128 *z = xF->data.v;
    __m128 *r = xF->data.v + sizeq;
    __m128 *hbar = xF->data.v + sizeq + sizeq;
    for (size_t i = 0; i < sizeq; i++) {
        r[i] *= istate->data.v[i];
    }
    cblas_sgemv(CblasColMajor, CblasTrans, sW2->nr, sW2->nc, 1.0, sW2->data.f,
                sW2->stride, (float *)r, 1, 1.0, (float *)hbar, 1);
    for (size_t i = 0; i < sizeq; i++) {
        hbar[i] = TANHFV(hbar[i]);
    }

    const __m128 ones = _mm_set1_ps(1.0f);
    for (size_t i = 0; i < sizeq ; i++) {
        ostate->data.v[i] = z[i] * istate->data.v[i] + (ones - z[i]) * hbar[i];
    }
}


flappie_matrix grumod_forward(const_flappie_matrix X, const_flappie_matrix sW,
                               flappie_matrix ostate) {
    RETURN_NULL_IF(NULL == X, NULL);

    assert(NULL != sW);

    const size_t bsize = X->nc;
    const size_t size = sW->nr;
    assert(X->nr == 3 * size);
    assert(sW->nc == 3 * size);

    ostate = remake_flappie_matrix(ostate, size, bsize);
    RETURN_NULL_IF(NULL == ostate, NULL);

    flappie_matrix tmp = make_flappie_matrix(3 * size, 1);
    if(NULL == tmp){
        //  Memory allocation falled, clean-up and return
        free(ostate);
        return NULL;
    }

    /* First step state is zero.  Set second column of ostate to zero and use that */
    _Mat xCol, sCol1, sCol2;
    memset(ostate->data.v + ostate->nrq, 0, ostate->nrq * sizeof(__m128));
    xCol = *X;
    sCol1 = *ostate;
    sCol2 = *ostate;
    xCol.nc = sCol1.nc = sCol2.nc = 1;
    sCol1.data.v = ostate->data.v + ostate->nrq;
    sCol2.data.v = ostate->data.v;
    grumod_step(&xCol, &sCol1, sW, tmp, &sCol2);
    for (size_t i = 1; i < bsize; i++) {
        xCol.data.v = X->data.v + i * X->nrq;
        sCol1.data.v = ostate->data.v + (i - 1) * ostate->nrq;
        sCol2.data.v = ostate->data.v + i * ostate->nrq;
        grumod_step(&xCol, &sCol1, sW, tmp, &sCol2);
    }

    tmp = free_flappie_matrix(tmp);

    assert(validate_flappie_matrix
           (ostate, -1.0, 1.0, 0.0, true, __FILE__, __LINE__));
    return ostate;
}


flappie_matrix grumod_backward(const_flappie_matrix X, const_flappie_matrix sW,
                                flappie_matrix ostate) {
    RETURN_NULL_IF(NULL == X, NULL);
    assert(NULL != sW);

    const size_t size = sW->nr;
    const size_t bsize = X->nc;
    assert(X->nr == 3 * size);
    assert(sW->nc == 3 * size);

    ostate = remake_flappie_matrix(ostate, size, bsize);
    RETURN_NULL_IF(NULL == ostate, NULL);

    flappie_matrix tmp = make_flappie_matrix(3 * size, 1);
    if(NULL == tmp){
        //  Memory allocation falled, clean-up and return
        free(ostate);
        return NULL;
    }

    /* First step state is zero.  Set first column of ostate to zero and use that */
    _Mat xCol, sCol1, sCol2;
    memset(ostate->data.v, 0, ostate->nrq * sizeof(__m128));
    xCol = *X;
    sCol1 = *ostate;
    sCol2 = *ostate;
    xCol.nc = sCol1.nc = sCol2.nc = 1;
    xCol.data.v = X->data.v + (X->nc - 1) * X->nrq;
    sCol1.data.v = ostate->data.v;
    sCol2.data.v = ostate->data.v + (ostate->nc - 1) * ostate->nrq;
    grumod_step(&xCol, &sCol1, sW, tmp, &sCol2);
    for (size_t i = 1; i < bsize; i++) {
        const size_t index = bsize - i - 1;
        xCol.data.v = X->data.v + index * X->nrq;
        sCol1.data.v = ostate->data.v + (index + 1) * ostate->nrq;
        sCol2.data.v = ostate->data.v + index * ostate->nrq;
        grumod_step(&xCol, &sCol1, sW, tmp, &sCol2);
    }

    tmp = free_flappie_matrix(tmp);

    assert(validate_flappie_matrix
           (ostate, -1.0, 1.0, 0.0, true, __FILE__, __LINE__));
    return ostate;
}


void grumod_step(const_flappie_matrix x, const_flappie_matrix istate,
                 const_flappie_matrix sW, flappie_matrix xF,
                 flappie_matrix ostate) {
    /* Perform a single modified GRU step
     * x      is [isize]
     * istate is [size]
     * xW     is [isize, 3 * size]
     * sW     is [size, 2 * size]
     * sW2    is [size, size]
     * bias   is [3 * size]
     * xF     is [3 * size]
     * ostate is [size]
     */
    assert(NULL != x);
    assert(NULL != sW);
    const size_t size = istate->nr;
    assert(x->nr == 3 * size);
    assert(size % 4 == 0);  // Vectorisation assumes size divisible by 4
    const size_t sizeq = size / 4;
    assert(size == sW->nr);
    assert(3 * size == sW->nc);
    assert(3 * size == xF->nr);
    assert(size == ostate->nr);


    // Copy input vector = iW x + b to temporary vector and zero last chunk
    memcpy(xF->data.v, x->data.v, x->nrq * sizeof(__m128));
    memset(xF->data.v + sizeq + sizeq, 0, sizeq *sizeof(__m128));
    /*  Add sW * istate to first 3 * size elts of xF
     *  then apply gate function to get r and z
     */
    cblas_sgemv(CblasColMajor, CblasTrans, sW->nr, sW->nc, 1.0, sW->data.f,
                sW->stride, istate->data.f, 1, 1.0, xF->data.f, 1);
    for (size_t i = 0; i < (sizeq + sizeq); i++) {
        xF->data.v[i] = LOGISTICFV(xF->data.v[i]);
    }

    const __m128 *z = xF->data.v;
    const __m128 *r = xF->data.v + sizeq;
    __m128 *hbar = xF->data.v + sizeq + sizeq;
    for (size_t i = 0; i < sizeq; i++) {
        hbar[i] = r[i] * hbar[i] + x->data.v[sizeq + sizeq + i];
    }
    for (size_t i = 0; i < sizeq; i++) {
        hbar[i] = TANHFV(hbar[i]);
    }

    const __m128 ones = _mm_set1_ps(1.0f);
    for (size_t i = 0; i < sizeq ; i++) {
        ostate->data.v[i] = z[i] * istate->data.v[i] + (ones - z[i]) * hbar[i];
    }
}


flappie_matrix gru_relu_forward(const_flappie_matrix X, const_flappie_matrix sW,
                                const_flappie_matrix sW2, flappie_matrix ostate) {
    RETURN_NULL_IF(NULL == X, NULL);

    assert(NULL != sW);
    assert(NULL != sW2);

    const size_t bsize = X->nc;
    const size_t size = sW2->nc;
    assert(X->nr == 3 * size);
    assert(sW->nr == size);
    assert(sW2->nr == size);
    assert(sW->nc == 2 * size);
    assert(sW2->nc == size);

    ostate = remake_flappie_matrix(ostate, size, bsize);
    RETURN_NULL_IF(NULL == ostate, NULL);

    flappie_matrix tmp = make_flappie_matrix(3 * size, 1);
    if(NULL == tmp){
        //  Memory allocation falled, clean-up and return
        free(ostate);
        return NULL;
    }

    /* First step state is zero.  Set second column of ostate to zero and use that */
    _Mat xCol, sCol1, sCol2;
    memset(ostate->data.v + ostate->nrq, 0, ostate->nrq * sizeof(__m128));
    xCol = *X;
    sCol1 = *ostate;
    sCol2 = *ostate;
    xCol.nc = sCol1.nc = sCol2.nc = 1;
    sCol1.data.v = ostate->data.v + ostate->nrq;
    sCol2.data.v = ostate->data.v;
    gru_relu_step(&xCol, &sCol1, sW, sW2, tmp, &sCol2);
    for (size_t i = 1; i < bsize; i++) {
        xCol.data.v = X->data.v + i * X->nrq;
        sCol1.data.v = ostate->data.v + (i - 1) * ostate->nrq;
        sCol2.data.v = ostate->data.v + i * ostate->nrq;
        gru_relu_step(&xCol, &sCol1, sW, sW2, tmp, &sCol2);
    }

    tmp = free_flappie_matrix(tmp);

    assert(validate_flappie_matrix
           (ostate, 0.0, HUGE_VAL, 0.0, true, __FILE__, __LINE__));
    return ostate;
}


flappie_matrix gru_relu_backward(const_flappie_matrix X, const_flappie_matrix sW,
                                 const_flappie_matrix sW2, flappie_matrix ostate) {
    RETURN_NULL_IF(NULL == X, NULL);
    assert(NULL != sW);
    assert(NULL != sW2);

    const size_t size = sW2->nc;
    const size_t bsize = X->nc;
    assert(X->nr == 3 * size);
    assert(sW->nr == size);
    assert(sW2->nr == size);
    assert(sW->nc == 2 * size);
    assert(sW2->nc == size);

    ostate = remake_flappie_matrix(ostate, size, bsize);
    RETURN_NULL_IF(NULL == ostate, NULL);

    flappie_matrix tmp = make_flappie_matrix(3 * size, 1);
    if(NULL == tmp){
        //  Memory allocation falled, clean-up and return
        free(ostate);
        return NULL;
    }

    /* First step state is zero.  Set first column of ostate to zero and use that */
    _Mat xCol, sCol1, sCol2;
    memset(ostate->data.v, 0, ostate->nrq * sizeof(__m128));
    xCol = *X;
    sCol1 = *ostate;
    sCol2 = *ostate;
    xCol.nc = sCol1.nc = sCol2.nc = 1;
    xCol.data.v = X->data.v + (X->nc - 1) * X->nrq;
    sCol1.data.v = ostate->data.v;
    sCol2.data.v = ostate->data.v + (ostate->nc - 1) * ostate->nrq;
    gru_relu_step(&xCol, &sCol1, sW, sW2, tmp, &sCol2);
    for (size_t i = 1; i < bsize; i++) {
        const size_t index = bsize - i - 1;
        xCol.data.v = X->data.v + index * X->nrq;
        sCol1.data.v = ostate->data.v + (index + 1) * ostate->nrq;
        sCol2.data.v = ostate->data.v + index * ostate->nrq;
        gru_relu_step(&xCol, &sCol1, sW, sW2, tmp, &sCol2);
    }

    tmp = free_flappie_matrix(tmp);

    assert(validate_flappie_matrix
           (ostate, 0.0, HUGE_VAL, 0.0, true, __FILE__, __LINE__));
    return ostate;
}


void gru_relu_step(const_flappie_matrix x, const_flappie_matrix istate,
                   const_flappie_matrix sW, const_flappie_matrix sW2,
                   flappie_matrix xF, flappie_matrix ostate) {
    /* Perform a single GRU step
     * x      is [isize]
     * istate is [size]
     * xW     is [isize, 3 * size]
     * sW     is [size, 2 * size]
     * sW2    is [size, size]
     * bias   is [3 * size]
     * xF     is [3 * size]
     * ostate is [size]
     */
    assert(NULL != x);
    assert(NULL != sW);
    assert(NULL != sW2);
    const size_t size = istate->nr;
    assert(x->nr == 3 * size);
    assert(size % 4 == 0);  // Vectorisation assumes size divisible by 4
    const size_t sizeq = size / 4;
    assert(size == sW->nr);
    assert(2 * size == sW->nc);
    assert(size == sW2->nr);
    assert(size == sW2->nc);
    assert(3 * size == xF->nr);
    assert(size == ostate->nr);


    // Copy input vector = iW x + b to temporary vector
    memcpy(xF->data.v, x->data.v, x->nrq * sizeof(__m128));
    /*  Add sW * istate to first 2 * size elts of xF
     *  then apply gate function to get r and z
     */
    cblas_sgemv(CblasColMajor, CblasTrans, sW->nr, sW->nc, 1.0, sW->data.f,
                sW->stride, istate->data.f, 1, 1.0, xF->data.f, 1);
    for (size_t i = 0; i < (sizeq +sizeq); i++) {
        xF->data.v[i] = LOGISTICFV(xF->data.v[i]);
    }

    const __m128 *z = xF->data.v;
    __m128 *r = xF->data.v + sizeq;
    __m128 *hbar = xF->data.v + sizeq + sizeq;
    for (size_t i = 0; i < sizeq; i++) {
        r[i] *= istate->data.v[i];
    }
    cblas_sgemv(CblasColMajor, CblasTrans, sW2->nr, sW2->nc, 1.0, sW2->data.f,
                sW2->stride, (float *)r, 1, 1.0, (float *)hbar, 1);
    for (size_t i = 0; i < sizeq; i++) {
        hbar[i] = relufv(hbar[i]);
    }

    const __m128 ones = _mm_set1_ps(1.0f);
    for (size_t i = 0; i < sizeq ; i++) {
        ostate->data.v[i] = z[i] * istate->data.v[i] + (ones - z[i]) * hbar[i];
    }
}


flappie_matrix lstm_forward(const_flappie_matrix Xaffine,
                             const_flappie_matrix sW,
                             flappie_matrix output) {
    RETURN_NULL_IF(NULL == Xaffine, NULL);
    assert(NULL != sW);

    const size_t size = sW->nr;
    const size_t bsize = Xaffine->nc;
    assert(Xaffine->nr == 4 * size);
    assert(sW->nc == 4 * size);

    output = remake_flappie_matrix(output, size, bsize);
    RETURN_NULL_IF(NULL == output, NULL);

    flappie_matrix tmp = make_flappie_matrix(4 * size, 1);
    flappie_matrix state = make_flappie_matrix(size, 1);
    if(NULL == tmp || NULL == state){
        //  Memory allocation falled, clean-up and return
        free(state);
        free(tmp);
        free(output);
        return NULL;
    }

    /* First step state & output are zero.  Set second column of output to zero and use that */
    memset(output->data.v + output->nrq, 0, output->nrq * sizeof(__m128));
    _Mat xCol, sCol1, sCol2;
    xCol = *Xaffine;
    sCol1 = *output;
    sCol2 = *output;
    xCol.nc = sCol1.nc = sCol2.nc = 1;
    sCol1.data.v = output->data.v + output->nrq;
    sCol2.data.v = output->data.v;
    lstm_step(&xCol, &sCol1, sW, tmp, state, &sCol2);
    for (size_t i = 1; i < bsize; i++) {
        xCol.data.v = Xaffine->data.v + i * Xaffine->nrq;
        sCol1.data.v = output->data.v + (i - 1) * output->nrq;
        sCol2.data.v = output->data.v + i * output->nrq;
        lstm_step(&xCol, &sCol1, sW, tmp, state, &sCol2);
    }

    state = free_flappie_matrix(state);
    tmp = free_flappie_matrix(tmp);

    assert(validate_flappie_matrix
           (output, -1.0, 1.0, 0.0, true, __FILE__, __LINE__));
    return output;
}


flappie_matrix lstm_backward(const_flappie_matrix Xaffine,
                              const_flappie_matrix sW,
                              flappie_matrix output) {
    RETURN_NULL_IF(NULL == Xaffine, NULL);
    assert(NULL != sW);

    const size_t size = sW->nr;
    const size_t bsize = Xaffine->nc;
    assert(Xaffine->nr == 4 * size);
    assert(sW->nc == 4 * size);

    output = remake_flappie_matrix(output, size, bsize);
    RETURN_NULL_IF(NULL == output, NULL);

    flappie_matrix tmp = make_flappie_matrix(4 * size, 1);
    flappie_matrix state = make_flappie_matrix(size, 1);
    if(NULL == tmp || NULL == state){
        //  Memory allocation falled, clean-up and return
        free(state);
        free(tmp);
        free(output);
        return NULL;
    }

    /* First step state is zero.  Set first column of ostate to zero and use that */
    memset(output->data.v, 0, output->nrq * sizeof(__m128));
    _Mat xCol, sCol1, sCol2;
    xCol = *Xaffine;
    sCol1 = *output;
    sCol2 = *output;
    xCol.nc = sCol1.nc = sCol2.nc = 1;
    xCol.data.v = Xaffine->data.v + (bsize - 1) * Xaffine->nrq;
    sCol1.data.v = output->data.v;
    sCol2.data.v = output->data.v + (bsize - 1) * output->nrq;
    lstm_step(&xCol, &sCol1, sW, tmp, state, &sCol2);
    for (size_t i = 1; i < bsize; i++) {
        const size_t index = bsize - i - 1;
        xCol.data.v = Xaffine->data.v + index * Xaffine->nrq;
        sCol1.data.v = output->data.v + (index + 1) * output->nrq;
        sCol2.data.v = output->data.v + index * output->nrq;
        lstm_step(&xCol, &sCol1, sW, tmp, state, &sCol2);
    }

    state = free_flappie_matrix(state);
    tmp = free_flappie_matrix(tmp);

    assert(validate_flappie_matrix
           (output, -1.0, 1.0, 0.0, true, __FILE__, __LINE__));
    return output;
}


void lstm_step(const_flappie_matrix xAffine, const_flappie_matrix out_prev,
               const_flappie_matrix sW,
               flappie_matrix xF, flappie_matrix state,
               flappie_matrix output) {
    /* Perform a single LSTM step
     * xAffine  is [isize] (== iW x + b, where x is the input to the LSTM layer)
     * out_prev is [size]
     * sW       is [size, 4 * size]
     * peep     is [4 * size]
     * xF       is [4 * size]
     * state    is [size]
     * output   is [size]
     */
    assert(NULL != xAffine);
    assert(NULL != out_prev);
    assert(NULL != sW);
    assert(NULL != xF);
    assert(NULL != state);
    assert(NULL != output);
    const size_t size = state->nr;
    assert(xAffine->nr == 4 * size);
    assert(size == out_prev->nr);
    assert(size == sW->nr);
    assert(4 * size == sW->nc);
    assert(4 * size == xF->nr);
    assert(size == output->nr);

    // Copy input vector = iW x + b to temporary vector
    memcpy(xF->data.v, xAffine->data.v, xAffine->nrq * sizeof(__m128));
    //  + sW' * xprev
    cblas_sgemv(CblasColMajor, CblasTrans, sW->nr, sW->nc, 1.0, sW->data.f,
                sW->stride, out_prev->data.f, 1, 1.0, xF->data.f, 1);

    assert(size % 4 == 0);  // Vectorisation assumes size divisible by 4
    const size_t sizeq = size / 4;
    for (size_t i = 0; i < sizeq; i++) {
        // Forget gate
        __m128 forget = LOGISTICFV(xF->data.v[sizeq + i])
                      * state->data.v[i];
        // Update gate
        __m128 update = LOGISTICFV(xF->data.v[i])
                      * TANHFV(xF->data.v[2 * sizeq + i]);
        state->data.v[i] = _mm_add_ps(forget, update);
        // Output gate
        output->data.v[i] = LOGISTICFV(xF->data.v[3 * sizeq + i])
                          * TANHFV(state->data.v[i]);
    }
}


size_t nbase_from_flipflop_nparam(size_t nparam){
    size_t nbase = roundf((-1.0f + sqrtf(1 + 2 * nparam)) / 2.0f);
    return nbase;
}


double crf_manystay_partition_function(const_flappie_matrix C){
    RETURN_NULL_IF(NULL == C, NAN);
    const size_t nbase = nbase_from_flipflop_nparam(C->nr);
    const size_t nstate = nbase + nbase;
    assert(nstate * (nbase + 1) == C->nr);

    double * mem = calloc(2 * nstate, sizeof(double));
    RETURN_NULL_IF(NULL==mem, NAN);

    double * curr = mem;
    double * prev = mem + nstate;

    for(size_t c=0 ; c < C->nc ; c++){
        const size_t offset = c * C->stride;
        const size_t offset_stay = offset + nstate * nbase;
        //  Swap
        {
            double * tmp = curr;
            curr = prev;
            prev = tmp;
        }

        for(size_t stay=nbase ; stay < nstate ; stay++){
           const size_t from_base = stay - nbase;
           curr[stay] = logsumexp(prev[stay] + C->data.f[offset_stay + stay],
                                   prev[from_base] + C->data.f[offset_stay + from_base]);
        }
        for(size_t to_state=0 ; to_state < nbase ; to_state++){
            const size_t offsetS = offset + to_state * nstate;
            curr[to_state] = C->data.f[offsetS + 0] + prev[0];
            for(size_t from_state=1 ; from_state < nstate ; from_state++){
                curr[to_state] = logsumexp(curr[to_state], C->data.f[offsetS + from_state] + prev[from_state]);
            }
        }
    }

    double logZ = curr[0];
    for(size_t st=1 ; st < nstate ; st++){
        logZ = logsumexp(logZ, curr[st]);
    }

    free(mem);

    return logZ;
}


flappie_matrix globalnorm_manystay(const_flappie_matrix X, const_flappie_matrix W,
                                    const_flappie_matrix b, float temperature, flappie_matrix C) {
    C = affine_map(X, W, b, C);
    RETURN_NULL_IF(NULL == C, NULL);
    tanh_activation_inplace(C);
    shift_scale_matrix_inplace(C, 0.0f, temperature / 5.0f);

    float logZ = crf_manystay_partition_function(C) / (double)C->nc;

    for(size_t c=0 ; c < C->nc ; c++){
        const size_t offset = c * C->stride;
        for(size_t r=0 ; r < C->nr ; r++){
            C->data.f[offset + r] -= logZ;
        }
    }


    return C;
}


flappie_matrix globalnorm_flipflop(const_flappie_matrix X, const_flappie_matrix W,
                                    const_flappie_matrix b, float temperature, flappie_matrix C) {
    return globalnorm_manystay(X, W, b, temperature, C);
}


/**  Calculates number of bases
 *
 *   @param nparams
 *
 *   @returns Number of bases
 **/
size_t nbase_from_runlength_nparam(size_t nparam){
    size_t nbase = nparam / 4;
    assert(4 * nbase == nparam);
    return nbase;
}

/**  Partition function for Run-length encoded model
 *
 *   @param C Transition parameters runlength model
 *
 *   @returns Logarithm of partition function
 **/
double runlength_partition_function(const_flappie_matrix C){
    RETURN_NULL_IF(NULL == C, NAN);

    const size_t nparam = C->nr;
    const size_t nbase = nbase_from_runlength_nparam(nparam);

    double * mem = calloc(2 * nbase, sizeof(double));
    RETURN_NULL_IF(NULL==mem, NAN);

    double * curr = mem;
    double * prev = mem + nbase;

    for(size_t c=0 ; c < C->nc ; c++){
        const size_t offset_move = c * C->stride + nbase + nbase;
        const size_t offset_stay = offset_move + nbase;
        //  Swap
        {
            double * tmp = curr;
            curr = prev;
            prev = tmp;
        }

        for(size_t b1=0 ; b1 < nbase ; b1++){
            curr[b1] = -HUGE_VAL;
            // Move from different base
            for(size_t b2=0 ; b2 < nbase ; b2++){
                if(b1 != b2){
                    curr[b1] = logsumexp(curr[b1], prev[b2]);
                }
            }
            curr[b1] += C->data.f[offset_move + b1];
        }
        for(size_t b=0 ; b < nbase ; b++){
            // Stay in same base
            curr[b] = logsumexp(curr[b], prev[b] + C->data.f[offset_stay + b]);
        }
    }


    double logZ = curr[0];
    for(size_t st=1 ; st < nbase ; st++){
        logZ = logsumexp(logZ, curr[st]);
    }

    free(mem);

    return logZ;
}


/**  Run-length encoded output layer
 *
 *   Performs initial linear transform and then scales all parameters appropriately.
 *
 *   Shape parameters
 *        x -> 1 + softplus(x)
 *   Scale parameters
 *        x -> ETA + softplus(x)
 *   Transition parameters
 *        x -> 5 tanh(x), global normalisation over all x
 *
 *
 *   @param X Input to layer
 *   @param W Weights for initial linear transform
 *   @param b Bias for initial linear transform
 *   @param temperature Temperature to normalise transition weights at
 *   @param C Flappie to write output parameters into.  Allocated if NULL.
 *
 *   @returns Parameters for run-length encoded model
 **/
flappie_matrix globalnorm_runlength(const_flappie_matrix X, const_flappie_matrix W,
                                    const_flappie_matrix b, float temperature,
                                    flappie_matrix C){
    const float ETA = 1e-1f;
    const size_t nbase = b->nr / 4;
    assert(nbase * 4 == b->nr);
    C = affine_map(X, W, b, C);
    RETURN_NULL_IF(NULL == C, NULL);

    for(size_t c=0 ; c < C->nc ; c++){
        const size_t offset = c * C->stride;
        for(size_t b=0 ; b < nbase ; b++){
            C->data.f[offset + b] = 1.0f + softplusf(C->data.f[offset + b]);
            C->data.f[offset + nbase + b] = ETA + softplusf(C->data.f[offset + nbase + b]);
            C->data.f[offset + 2 * nbase + b] = 5.0f * tanhf(C->data.f[offset + 2 * nbase + b]) / temperature;
            C->data.f[offset + 3 * nbase + b] = 5.0f * tanhf(C->data.f[offset + 3 * nbase + b]) / temperature;
        }
    }

    float logZ = runlength_partition_function(C) / (float)C->nc;

    for(size_t c=0 ; c < C->nc ; c++){
        const size_t offset = c * C->stride + 2 * nbase;
        for(size_t r=0 ; r < 2 * nbase ; r++){
            C->data.f[offset + r] -= logZ;
        }
    }

    return C;
}


/**  Calculates number of bases
 *
 *   @param nparams
 *
 *   @returns Number of bases
 **/
size_t nbase_from_crf_runlength_nparam(size_t nparam){
    // By numerical accident, answer is same as flipflip
    // 2 * nbase * nbase + 2 * nbase
    return nbase_from_flipflop_nparam(nparam);
}

static inline size_t rle_trans_lookup(size_t base_from, bool stay_from,
                                      size_t base_to, bool stay_to,
                                      size_t nbase){
    assert(stay_to ^ (base_from != base_to));
    return base_to * 2 * nbase + base_from + (stay_from ? nbase : 0);
}


/**  Partition function for new version Run-length encoded model
 *
 *   @param C Transition parameters runlength model
 *
 *   @returns Logarithm of partition function
 **/
double runlengthV2_partition_function(const_flappie_matrix C){
    RETURN_NULL_IF(NULL == C, NAN);

    const size_t nparam = C->nr;
    const size_t nbase = nbase_from_crf_runlength_nparam(nparam);
    const size_t nstate = nbase + nbase;

    double * mem = calloc(2 * nstate, sizeof(double));
    RETURN_NULL_IF(NULL==mem, NAN);

    double * curr = mem;
    double * prev = mem + nstate;

    for(size_t c=0 ; c < C->nc ; c++){
        const size_t offset = c * C->stride + nstate;
        //  Swap
        {
            double * tmp = curr;
            curr = prev;
            prev = tmp;
        }

        for(size_t b1=0 ; b1 < nbase ; b1++){
            curr[b1] = -HUGE_VAL;
            for(size_t b2=0 ; b2 < nbase ; b2++){
                if(b1 == b2){
                    continue;
                }
                // Move from different base (move)
                curr[b1] = logsumexp(curr[b1], prev[b2] + C->data.f[offset + rle_trans_lookup(b2, false, b1, false, nbase)]);
                // Move from different base (stay)
                curr[b1] = logsumexp(curr[b1], prev[b2 + nbase] + C->data.f[offset + rle_trans_lookup(b2, true, b1, false, nbase)]);
            }
            for(size_t b=0 ; b < nbase ; b++){
                // Stay in same base
                curr[b + nbase] = logsumexpf(prev[b] + C->data.f[offset + rle_trans_lookup(b, false, b, true, nbase)],
                                             prev[b + nbase] + C->data.f[offset + rle_trans_lookup(b, true, b, true, nbase)]);
            }
        }
    }

    double logZ = curr[0];
    for(size_t st=1 ; st < nstate ; st++){
        logZ = logsumexp(logZ, curr[st]);
    }

    free(mem);

    return logZ;
}


/**  Run-length encoded output layer (new version)
 *
 *   Performs initial linear transform and then scales all parameters appropriately.
 *
 *   Shape parameters (4 parameters per block)
 *        x -> 1 + softplus(x)
 *   Scale parameters (4 parameters pe block)
 *        x -> ETA + softplus(x)
 *   Transition parameters (32 parameters per block)
 *        x -> 5 tanh(x), global normalisation over all x
 *
 *
 *   @param X Input to layer
 *   @param W Weights for initial linear transform
 *   @param b Bias for initial linear transform
 *   @param temperature Temperature to normalise transition weights at
 *   @param C Flappie to write output parameters into.  Allocated if NULL.
 *
 *   @returns Parameters for run-length encoded model
 **/
flappie_matrix globalnorm_runlengthV2(const_flappie_matrix X, const_flappie_matrix W,
                                      const_flappie_matrix b, float temperature,
                                      flappie_matrix C){
    C = affine_map(X, W, b, C);
    RETURN_NULL_IF(NULL == C, NULL);
    const size_t nbase = nbase_from_crf_runlength_nparam(C->nr);
    const size_t nrunparam = nbase + nbase;

    for(size_t c=0 ; c < C->nc ; c++){
        const size_t offset = c * C->stride;
        for(size_t b=0 ; b < nbase ; b++){
            //  Shift and scale parameters
            C->data.f[offset + b] = 1.0f + softplusf(C->data.f[offset + b]);
            C->data.f[offset + nbase + b] = 1e-8f + softplusf(C->data.f[offset + nbase + b]);
        }
        for(size_t param=nrunparam ; param < C->nr ; param++){
            //  Transition parameters
            C->data.f[offset + param] = 5.0f * tanhf(C->data.f[offset + param]) / temperature;
        }
    }

    const float logZ = runlengthV2_partition_function(C) / (float)C->nc;

    for(size_t c=0 ; c < C->nc ; c++){
        const size_t offset = c * C->stride;
        for(size_t r=nrunparam ; r < C->nr ; r++){
            C->data.f[offset + r] -= logZ;
        }
    }

    return C;
}


