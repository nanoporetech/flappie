/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#define BANANA 1
#include <CUnit/CUnit.h>
#include <stdbool.h>
#include <stdlib.h>

#include <layers.h>
#include "test_common.h"

static float test_conv_tol = 1e-5;

typedef struct {
    float *elt;
    size_t len;
} Vec;

static float _xrange_odd[] =
    { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
static Vec const xrange_odd = {.elt = _xrange_odd,.len = 11 };

static float _xrange_even[] =
    { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
static Vec const xrange_even = {.elt = _xrange_even,.len = 10 };

bool compare_vecs_for_equality(Vec const x, Vec const y, float tol) {
    if (x.len != y.len) {
        return false;
    }
    for (size_t i = 0; i < x.len; ++i) {
        if (fabsf(x.elt[i] - y.elt[i]) > tol) {
            return false;
        }
    }
    return true;
}

/**  Simple implementation of convolution
 *
 *   @param x Vector to form convolution over.
 *   @param f Filter for convolution
 *   @param winlen Window of convolution( (length of filter)
 *   @returns vector containing filtered data with same length as x
 **/
Vec simple_convolution(Vec const x, Vec const f) {
    const size_t winlen = f.len;

    CU_ASSERT_PTR_NOT_NULL_FATAL(x.elt);
    CU_ASSERT(x.len > 0);
    CU_ASSERT_PTR_NOT_NULL_FATAL(f.elt);
    CU_ASSERT(winlen > 0);
    CU_ASSERT(x.len >= winlen);

    // Pad input.
    const size_t padL = (winlen - 1) / 2;
    const size_t padR = winlen / 2;
    float *xpad = calloc(padL + padR + x.len, sizeof(float));
    CU_ASSERT_PTR_NOT_NULL_FATAL(xpad);
    memcpy(xpad + padL, x.elt, x.len * sizeof(float));

    // Vector for output
    Vec y = { NULL, x.len };
    y.elt = calloc(x.len, sizeof(float));
    CU_ASSERT_PTR_NOT_NULL_FATAL(y.elt);

    for (size_t start = 0; start < x.len; ++start) {
        for (size_t w = 0; w < winlen; ++w) {
            y.elt[start] += f.elt[w] * xpad[start + w];
        }
    }

    free(xpad);
    return y;
}


/**  Created strided vector from full vector
 *
 *   @param v Vector from which to create strided vector
 *   @param stride Length of stride
 *
 *   @returns vector containing stride
 *
 **/
Vec simple_stride(Vec const v, int stride) {
    CU_ASSERT_PTR_NOT_NULL_FATAL(v.elt);
    CU_ASSERT(v.len > 0);
    CU_ASSERT(stride > 0);

    const size_t newlen = (v.len + stride - 1) / stride;
    Vec y = { NULL, newlen };
    y.elt = calloc(newlen, sizeof(float));
    CU_ASSERT_PTR_NOT_NULL_FATAL(y.elt);

    for (size_t i = 0, j = 0; i < v.len; i += stride, ++j) {
        y.elt[j] = v.elt[i];
    }

    return y;
}


/**  Initialise test
 *
 *   @returns 0 on success, non-zero on failure
 **/
static flappie_matrix mat_even = NULL;
static flappie_matrix mat_odd = NULL;
int init_test_convolution(void) {
    mat_odd = mat_from_array(_xrange_odd, 1, 11);
    if (NULL == mat_odd) {
        return 1;
    }
    mat_even = mat_from_array(_xrange_even, 1, 10);
    if (NULL == mat_even) {
        mat_odd = free_flappie_matrix(mat_odd);
        return 1;
    }
    return 0;
}


/**  Clean up after test
 *
 *   @returns 0 on success, non-zero on failure
 **/
int clean_test_convolution(void) {
    mat_even = free_flappie_matrix(mat_even);
    mat_odd = free_flappie_matrix(mat_odd);
    return 0;
}


/**
 *
 *   Firstly test that our ground-truth routines above give correct answers
 *
 **/

/**  Helper function to check application of stride to a given input/
 *
 *   @param input A vector to stride across
 *   @param expected A vector to compare result
 *   @param stride to apply
 **/
void check_stride_for_equality(Vec const input, Vec const expected, int stride) {
    Vec y = simple_stride(input, stride);
    CU_ASSERT_TRUE(compare_vecs_for_equality(y, expected, 0.0));
    free(y.elt);
}

void test_stride1_convolution(void) {
    check_stride_for_equality(xrange_odd, xrange_odd, 1);
    check_stride_for_equality(xrange_even, xrange_even, 1);
}

void test_stride2_convolution(void) {
    int const stride = 2;
    float _expected_odd[] = { 0.0, 2.0, 4.0, 6.0, 8.0, 10.0 };
    Vec const expected_odd = {
        .elt = _expected_odd,
        .len = 6
    };
    float _expected_even[] = { 0.0, 2.0, 4.0, 6.0, 8.0 };
    Vec const expected_even = {
        .elt = _expected_even,
        .len = 5
    };

    check_stride_for_equality(xrange_odd, expected_odd, stride);
    check_stride_for_equality(xrange_even, expected_even, stride);
}

void test_stride3_convolution(void) {
    int const stride = 3;
    float _expected_odd[] = { 0.0, 3.0, 6.0, 9.0 };
    Vec const expected_odd = {
        .elt = _expected_odd,
        .len = 4
    };
    float _expected_even[] = { 0.0, 3.0, 6.0, 9.0 };
    Vec const expected_even = {
        .elt = _expected_even,
        .len = 4
    };

    check_stride_for_equality(xrange_odd, expected_odd, stride);
    check_stride_for_equality(xrange_even, expected_even, stride);
}

void test_stride4_convolution(void) {
    int const stride = 4;
    float _expected_odd[] = { 0.0, 4.0, 8.0 };
    Vec const expected_odd = {
        .elt = _expected_odd,
        .len = 3
    };
    float _expected_even[] = { 0.0, 4.0, 8.0 };
    Vec const expected_even = {
        .elt = _expected_even,
        .len = 3
    };

    check_stride_for_equality(xrange_odd, expected_odd, stride);
    check_stride_for_equality(xrange_even, expected_even, stride);
}

void test_stride5_convolution(void) {
    int const stride = 5;
    float _expected_odd[] = { 0.0, 5.0, 10.0 };
    Vec const expected_odd = {
        .elt = _expected_odd,
        .len = 3
    };
    float _expected_even[] = { 0.0, 5.0 };
    Vec const expected_even = {
        .elt = _expected_even,
        .len = 2
    };

    check_stride_for_equality(xrange_odd, expected_odd, stride);
    check_stride_for_equality(xrange_even, expected_even, stride);
}

void check_convolution_for_equality(Vec const input, Vec const filter,
                                    Vec const expected, float tol) {
    Vec res = simple_convolution(input, filter);
    CU_ASSERT_PTR_NOT_NULL(res.elt);
    CU_ASSERT_TRUE(compare_vecs_for_equality(res, expected, tol));
    free(res.elt);
}

void test_convolution_ones_f1(void) {
    float _filter[] = { 1.0 };
    Vec const filter = {
        .elt = _filter,
        .len = 1
    };

    check_convolution_for_equality(xrange_odd, filter, xrange_odd,
                                   test_conv_tol);
    check_convolution_for_equality(xrange_even, filter, xrange_even,
                                   test_conv_tol);
}

void test_convolution_ones_f2(void) {
    float _filter[] = { 1.0, 1.0 };
    Vec const filter = {
        .elt = _filter,
        .len = 2
    };
    float _expected_odd[] =
        { 1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0, 10.0 };
    Vec const expected_odd = {
        .elt = _expected_odd,
        .len = 11
    };
    float _expected_even[] =
        { 1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 9.0 };
    Vec const expected_even = {
        .elt = _expected_even,
        .len = 10
    };

    check_convolution_for_equality(xrange_odd, filter, expected_odd,
                                   test_conv_tol);
    check_convolution_for_equality(xrange_even, filter, expected_even,
                                   test_conv_tol);
}

void test_convolution_ones_f3(void) {
    float _filter[] = { 1.0, 1.0, 1.0 };
    Vec const filter = {
        .elt = _filter,
        .len = 3
    };
    float _expected_odd[] =
        { 1.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 19.0 };
    Vec const expected_odd = {
        .elt = _expected_odd,
        .len = 11
    };
    float _expected_even[] =
        { 1.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 17.0 };
    Vec const expected_even = {
        .elt = _expected_even,
        .len = 10
    };

    check_convolution_for_equality(xrange_odd, filter, expected_odd,
                                   test_conv_tol);
    check_convolution_for_equality(xrange_even, filter, expected_even,
                                   test_conv_tol);
}

void test_convolution_ones_f4(void) {
    float _filter[] = { 1.0, 1.0, 1.0, 1.0 };
    Vec const filter = {
        .elt = _filter,
        .len = 4
    };
    float _expected_odd[] =
        { 3.0, 6.0, 10.0, 14.0, 18.0, 22.0, 26.0, 30.0, 34.0, 27.0, 19.0 };
    Vec const expected_odd = {
        .elt = _expected_odd,
        .len = 11
    };
    float _expected_even[] =
        { 3.0, 6.0, 10.0, 14.0, 18.0, 22.0, 26.0, 30.0, 24.0, 17.0 };
    Vec const expected_even = {
        .elt = _expected_even,
        .len = 10
    };

    check_convolution_for_equality(xrange_odd, filter, expected_odd,
                                   test_conv_tol);
    check_convolution_for_equality(xrange_even, filter, expected_even,
                                   test_conv_tol);
}

void test_convolution_ones_f5(void) {
    float _filter[] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
    Vec const filter = {
        .elt = _filter,
        .len = 5
    };
    float _expected_odd[] =
        { 3.0, 6.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 34.0, 27.0 };
    Vec const expected_odd = {
        .elt = _expected_odd,
        .len = 11
    };
    float _expected_even[] =
        { 3.0, 6.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 30.0, 24.0 };
    Vec const expected_even = {
        .elt = _expected_even,
        .len = 10
    };

    check_convolution_for_equality(xrange_odd, filter, expected_odd,
                                   test_conv_tol);
    check_convolution_for_equality(xrange_even, filter, expected_even,
                                   test_conv_tol);
}

void test_convolution_antisymmetric_f3(void) {
    float _filter[] = { -1.0, 0.0, 1.0 };
    Vec const filter = {
        .elt = _filter,
        .len = 3
    };
    float _expected_odd[] =
        { 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, -9.0 };
    Vec const expected_odd = {
        .elt = _expected_odd,
        .len = 11
    };
    float _expected_even[] =
        { 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, -8.0 };
    Vec const expected_even = {
        .elt = _expected_even,
        .len = 10
    };

    check_convolution_for_equality(xrange_odd, filter, expected_odd,
                                   test_conv_tol);
    check_convolution_for_equality(xrange_even, filter, expected_even,
                                   test_conv_tol);
}

void compare_convolution_to_baseline(const_flappie_matrix input,
                                     const_flappie_matrix filter,
                                     const_flappie_matrix bias, Vec input_base,
                                     Vec filter_base) {
    flappie_matrix res = convolution(input, filter, bias, 1, NULL);
    Vec res_base = simple_convolution(input_base, filter_base);
    CU_ASSERT_PTR_NOT_NULL_FATAL(res);
    CU_ASSERT_PTR_NOT_NULL_FATAL(res_base.elt);

    for (size_t i = 0; i < input_base.len; ++i) {
        CU_ASSERT_DOUBLE_EQUAL(res_base.elt[i], res->data.f[i * 4],
                               test_conv_tol);
    }

    free(res_base.elt);
    res = free_flappie_matrix(res);
}

void test_flappie_convolution_f1s1(void) {
    float _filter[] = { 1.0, 0.0, 0.0, 0.0 };
    float _bias[4] = { 0.0 };
    _Mat filter = {
        .nr = 1, .nrq = 1, .nc = 1, .stride=4,
        .data.f = _filter
    };
    _Mat bias = {
        .nr = 1, .nrq = 1, .nc = 1, .stride=4,
        .data.f = _bias
    };

    Vec filter_base = {
        .elt = _filter,
        .len = 1
    };

    compare_convolution_to_baseline(mat_odd, &filter, &bias, xrange_odd,
                                    filter_base);
    compare_convolution_to_baseline(mat_even, &filter, &bias, xrange_even,
                                    filter_base);
}


static const test_with_description tests[] = {
    {"Simple stride 1", test_stride1_convolution},
    {"Simple stride 2", test_stride2_convolution},
    {"Simple stride 3", test_stride3_convolution},
    {"Simple stride 4", test_stride4_convolution},
    {"Simple stride 5", test_stride5_convolution},
    {"Simple convolution, unit filter length 1", test_convolution_ones_f1},
    {"Simple convolution, unit filter length 2", test_convolution_ones_f2},
    {"Simple stride 1", test_stride1_convolution},
    {"Simple stride 2", test_stride2_convolution},
    {"Simple stride 3", test_stride3_convolution},
    {"Simple stride 4", test_stride4_convolution},
    {"Simple stride 5", test_stride5_convolution},
    {"Simple convolution, unit filter length 1", test_convolution_ones_f1},
    {"Simple convolution, unit filter length 2", test_convolution_ones_f2},
    {"Simple convolution, unit filter length 3", test_convolution_ones_f3},
    {"Simple convolution, unit filter length 4", test_convolution_ones_f4},
    {"Simple convolution, unit filter length 5", test_convolution_ones_f5},
    {"Simple convolution, antisymmetric filter length 3", test_convolution_antisymmetric_f3},
    {"Scrappie convolution, antisymmetric filter length 3", test_flappie_convolution_f1s1},
    {0}};

/**   Register tests with CUnit
 *
 *    @returns 0 on success, non-zero on failure
 **/
int register_test_convolution(void) {
    return flappie_register_test_suite("Convolution layer", init_test_convolution, clean_test_convolution, tests);
}
