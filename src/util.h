/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#pragma once
#ifndef UTIL_H
#    define UTIL_H
#    include <assert.h>
#    include <immintrin.h>
#    include <math.h>
#    include <stdbool.h>
#    include <stdint.h>
#    include "sse_mathfun.h"

#    ifdef FAST_LOG
#        define LOGFV fast_logfv
#    else
#        define LOGFV logfv
#    endif

#    ifdef FAST_EXP
#        define EXPFV fast_expfv
#    else
#        define EXPFV expfv
#    endif

#    ifdef FAST_TANH
#        define TANHFV fast_tanhfv
#    else
#        define TANHFV tanhfv
#    endif

#    ifdef FAST_LOGISTIC
#        define LOGISTICFV fast_logisticfv
#    else
#        define LOGISTICFV logisticfv
#    endif

#    ifdef FAST_ELU
#        define ELUFV fast_elufv
#    else
#        define ELUFV elufv
#    endif


/* From math.h */
#    ifndef M_LN2
#        define M_LN2          0.69314718055994530942  /* log_e 2 */
#    endif
#    ifndef M_LOG10E
#        define M_LOG10E       0.43429448190325182765  /* log_10 e */
#    endif

/* Create a vector of  ones.  */
extern __inline __m128 __attribute__ ((__gnu_inline__, __always_inline__))
    _mm_setone_ps(void) {
    return __extension__(__m128) {
    1.0f, 1.0f, 1.0f, 1.0f};
}

extern __inline __m128d __attribute__ ((__gnu_inline__, __always_inline__))
    _mm_setone_pd(void) {
    return __extension__(__m128d) {
    1.0, 1.0};
}


/**
 *   Standard implementations
 **/
static inline float logisticf(float x){
    return 1.0f / (1.0f + expf(-x));
}

static inline float eluf(float x){
    return (x >= 0.0f) ? x : expm1f(x);
}

static inline float softplusf(float x){
    return log1pf(expf(-fabsf(x))) + ((x >= 0.0f) ? x : 0.f);
}

static inline float powm1f(float x, float y){
    return expm1f(y * logf(x));
}


/**
 *   Laplace distribution and derivatives
 **/
static inline float loglaplace(float x, float loc, float sc, float logsc){
    return -fabsf(x - loc) / sc - logsc - M_LN2;
}

static inline float laplace(float x, float loc, float sc, float logsc){
    return expf(loglaplace(x, loc, sc, logsc));
}

static inline float dloglaplace_loc(float x, float loc, float sc){
    return ((x > loc)  - (x < loc)) / sc;
}

static inline float dloglaplace_scale(float x, float loc, float sc){
    return (fabsf(x - loc) / sc - 1.0) / sc;
}

static inline float dloglaplace_logscale(float x, float loc, float sc){
    return fabsf(x - loc) / sc - 1.0;
}

static inline float dlaplace_loc(float x, float loc, float sc, float logsc){
    return laplace(x, loc, sc, logsc) * dloglaplace_loc(x, loc, sc);
}

static inline float dlaplace_scale(float x, float loc, float sc, float logsc){
    return laplace(x, loc, sc, logsc) * dloglaplace_scale(x, loc, sc);
}

static inline float dlaplace_logscale(float x, float loc, float sc, float logsc){
    return laplace(x, loc, sc, logsc) * dloglaplace_logscale(x, loc, sc);
}

static inline float lchoosef(float n, float k){
        return lgammaf(n + 1.0f) - lgammaf(n - k + 1.0f) - lgammaf(k + 1.0f);
}

static inline float logdnegbinomf(float k, float r, float p){
    return k * logf(p) + r * log1pf(-p) + lchoosef(k + r - 1.0f, k);
}


/**
 *   Logistic distribution
 **/
static inline float plogisticf(float x){
    return 0.5f * (1.0f + tanhf(x / 2.0f));
}

static inline float logplogisticf(float x){
    return -log1pf(expf(-x));
}

static inline float qlogisticf(float p){
    return 2.0f * atanhf(2.0f * p - 1.0f);
}

static inline float dlogisticf(float x){
    const float p = plogisticf(x);
    return p * (1.0f - p);
}


/**
 *    Weibull distribution
 **/
static inline float pweibullf(float x, float sh, float sc){
    //  Cumulative density function of Weibull distribution
    assert(x >= 0.0);
    assert(sh > 0.0);
    assert(sc > 0.0);
    const float xsc = x / sc;
    const float p1 = powf(xsc, sh);
    return -expm1(-p1);
}

static inline float logpweibullf(float x, float sh, float sc){
    //  Log cumulative density function of Weibull distribution
    assert(x >= 0.0);
    assert(sh > 0.0);
    assert(sc > 0.0);
    const float xsc = x / sc;
    const float p1 = powf(xsc, sh);
    return logf(-expm1(-p1));
}

static inline float logcpweibullf(float x, float sh, float sc){
    //  Log complementary CDF of Weibull distribution
    assert(x >= 0.0);
    assert(sh > 0.0);
    assert(sc > 0.0);
    const float xsc = x / sc;
    const float p1 = powf(xsc, sh);
    return -p1;
}

static inline float dweibullf(float x, float sh, float sc){
    //  density of Weibull distribution
    assert(x >= 0.0);
    assert(sh > 0.0);
    assert(sc > 0.0);
    const float xsc = x / sc;
    const float p1 = powf(xsc, sh);
    return sh * p1 * expf(-p1) / x;
}


/**
 *   Discrete Weibull distribution
 **/
static inline float pdiscreteweibullf(float x, float sh, float sc){
    //  CDF of the Discrete Weibull distribution
    assert(x >= 0.0);
    assert(sh > 0.0);
    assert(sc > 0.0);
    return pweibullf(x + 1.0f, sh, sc);
}

static inline float logpdiscreteweibullf(float x, float sh, float sc){
    //  Log CDF of the Discrete Weibull distribution
    assert(x >= 0.0);
    assert(sh > 0.0);
    assert(sc > 0.0);
    return logpweibullf(x + 1.0f, sh, sc);
}

static inline float logcpdiscreteweibullf(float x, float sh, float sc){
    //  Log complementary CDF of the Discrete Weibull distribution
    assert(x >= 0.0);
    assert(sh > 0.0);
    assert(sc > 0.0);
    return logcpweibullf(x + 1.0f, sh, sc);
}

static inline float ddiscreteweibullf(float x, float sh, float sc){
    //  probability mass function of Discrete Weibull distribution
    assert(x >= 0.0);
    assert(sh > 0.0);
    assert(sc > 0.0);
    const float log_cprob1 = -powf(x / sc, sh);
    const float log_cprob2 = -powf((x + 1.0f) / sc, sh);
    const float delta_log_cprob = -log_cprob2 * powm1f(x / (1.0f + x), sh);
    return -expf(log_cprob1) * expm1f(delta_log_cprob);
}




// Constants for fast exp approximation.  See Schraudolph (1999)
#    define _A 12102203.161561485f
//  Minimum maximum relative error
//#    define _B 1064986822.5027076f
//  No bias at zero
#    define _B 1065353216.0f
//  Minimum RMS relative error
//#    define _B 1064866803.6193156f
//  Minimum mean relative error
//#    define _B 1064807268.0425934f
#    define _BOUND 88.02969193111305
static inline float fast_expf(float x) {
    x = fmaxf(-_BOUND, fminf(_BOUND, x));
    union {
        uint32_t i;
        float f;
    } value = {
    .i = (uint32_t) (_A * x + _B)};
    return value.f;
}

static inline float fast_logisticf(float x) {
    return 1.0 / (1.0 + fast_expf(-x));
}

static inline float fast_tanhf(float x) {
    const float y = fast_logisticf(x + x);
    return y + y - 1.0;
}

static inline float fast_eluf(float x){
    return (x >= 0.0f) ? x : (fast_expf(x) - 1.0);
}

static inline float logsumexpf(float x, float y){
        return fmaxf(x, y) + log1pf(expf(-fabsf(x-y)));
}

static inline double logsumexp(double x, double y){
        return fmax(x, y) + log1p(exp(-fabs(x-y)));
}


#define MAX_POST_PROB 0.99999
static inline float qscoref(float p){
    assert(p >= 0.0f && p <=1.0f);
    const float p_clip = (p < MAX_POST_PROB) ? p : MAX_POST_PROB;
    return -(10.0f * M_LOG10E) * log1pf(-p_clip);
}


static inline double qscore(double p){
    assert(p >= 0.0 && p <=1.0);
    const float p_clip = (p < MAX_POST_PROB) ? p : MAX_POST_PROB;
    return -(10.0 * M_LOG10E) * log1p(-p_clip);
}


static inline char phredf(float p){
    char ph = roundf(33.0f + qscoref(p));
    char ph_clip = (ph < 126) ? ph : 126;
    assert(ph_clip > 32 && ph_clip < 127);
    return ph_clip;
}


static inline char phred(float p){
    char ph = round(33.0 + qscore(p));
    char ph_clip = (ph < 126) ? ph : 126;
    assert(ph_clip > 32 && ph_clip < 127);
    return ph_clip;
}


/**
 *   Vectorised functions
 **/
static inline __m128 __attribute__ ((__always_inline__)) expfv(__m128 x) {
    __v4sf y = (__v4sf) x;
    return (__m128) exp_ps(y);
}

static inline __m128 __attribute__ ((__always_inline__)) logfv(__m128 x) {
    __v4sf y = (__v4sf) x;
    return (__m128) log_ps(y);
}

static inline __m128 __attribute__ ((__always_inline__)) logisticfv(__m128 x) {
    const __m128 ones = _mm_setone_ps();
    return _mm_div_ps(ones, _mm_add_ps(ones, expfv(-x)));
}

static inline __m128 __attribute__ ((__always_inline__)) tanhfv(__m128 x) {
    const __m128 y = logisticfv(_mm_add_ps(x, x));
    return _mm_sub_ps(_mm_add_ps(y, y), _mm_setone_ps());
}

static inline __m128 __attribute__ ((__always_inline__)) elufv(__m128 x) {
    if(0 == _mm_movemask_ps(x)){
        // All positive, early return.
        return x;
    }
    const __m128 mask = _mm_cmpge_ps(x, _mm_setzero_ps());
    const __m128 y = expfv(x) - _mm_setone_ps();
    return _mm_or_ps(_mm_and_ps(mask, x),  _mm_andnot_ps(mask, y));
}

static inline __m128 __attribute__ ((__always_inline__)) relufv(__m128 x) {
    if(0 == _mm_movemask_ps(x)){
        // All positive, early return.
        return x;
    }
    const __m128 mask = _mm_cmpge_ps(x, _mm_setzero_ps());
    return _mm_or_ps(_mm_and_ps(mask, x),  _mm_andnot_ps(mask, _mm_setzero_ps()));
}

/**
 *    Fast vectorised approximations
 **/
static inline __m128 fast_expfv(__m128 x) {
    const __m128 a = (__m128) (__v4sf) { _A, _A, _A, _A };
    const __m128 b = (__m128) (__v4sf) { _B, _B, _B, _B };
    const __m128 _bound = (__m128) (__v4sf) { _BOUND, _BOUND, _BOUND, _BOUND };
    x = _mm_max_ps(-_bound, _mm_min_ps(_bound, x));

    __m128 y = a * x + b;
    return _mm_castsi128_ps(_mm_cvtps_epi32(y));
}

static inline __m128 fast_logfv(__m128 x) {
#    define _Alogfv 8.262958294867817e-08f
#    define _Blogfv 1064872507.1541044f
    const __m128 a = (__m128) (__v4sf) { _Alogfv, _Alogfv, _Alogfv, _Alogfv };
    const __m128 b = (__m128) (__v4sf) { _Blogfv, _Blogfv, _Blogfv, _Blogfv };
    x = _mm_cvtepi32_ps(_mm_castps_si128(x));
    return a * (x - b);
}

static inline __m128 __attribute__ ((__always_inline__)) fast_logisticfv(__m128 x) {
    return _mm_rcp_ps(_mm_add_ps(_mm_setone_ps(), fast_expfv(-x)));
}

static inline __m128 __attribute__ ((__always_inline__)) fast_tanhfv(__m128 x) {
    const __m128 y = fast_logisticfv(_mm_add_ps(x, x));
    return _mm_sub_ps(_mm_add_ps(y, y), _mm_setone_ps());
}

static inline __m128 __attribute__ ((__always_inline__)) fast_elufv(__m128 x) {
    if(0 == _mm_movemask_ps(x)){
        // All positive, early return.
        return x;
    }
    const __m128 mask = _mm_cmpge_ps(x, _mm_setzero_ps());
    const __m128 y = fast_expfv(x) - _mm_setone_ps();
    return _mm_or_ps(_mm_and_ps(mask, x),  _mm_andnot_ps(mask, y));
}

int argmaxf(const float *x, size_t n);
int argminf(const float *x, size_t n);
float valmaxf(const float *x, size_t n);
float valminf(const float *x, size_t n);

static inline int iceil(int x, int y) {
    return (x + y - 1) / y;
}

static inline int ifloor(int x, int y) {
    return x / y;
}

static inline int imin(int x, int y){
    return (x < y) ? x : y;
}

static inline int imax(int x, int y){
    return (x > y) ? x : y;
}

void quantilef(const float *x, size_t nx, float *p, size_t np);
float medianf(const float *x, size_t n);
float madf(const float *x, size_t n, const float *med);
void medmad_normalise_array(float *x, size_t n);
void studentise_array_kahan(float *x, size_t n);
void difference_array(float *x, size_t n);
void filter_array(float *x, size_t n, float fill_val, float thresh);
void clip_array(float *x, size_t n, float thresh);
void qscore_array(float *x, size_t n);

bool equality_array(double const * x, double const * y, size_t n, double const tol);
bool equality_arrayf(float const * x, float const * y, size_t n, float const tol);
bool equality_arrayi(int const * x, int const * y, size_t n);
size_t ndiff_array(const int * x, size_t n);
bool mod_array_inplace(int * x, size_t n, int m);

char * phredstr_from_qscore(float * x, size_t n);

#endif                          /* UTIL_H */
