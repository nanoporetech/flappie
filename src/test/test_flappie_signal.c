/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#include <CUnit/Basic.h>
#include <err.h>
#include <stdbool.h>

#include "layers.h"
#include "flappie_common.h"
#include "flappie_structures.h"
#include "flappie_util.h"
#include "test_common.h"
#include "util.h"

static const char rawsignalfile[] = "raw_signal.crp";
static const char trimsignalfile[] = "trimmed_signal.crp";
static const char normsignalfile[] = "normalised_signal.crp";


static flappie_matrix rawsignal = NULL;
static flappie_matrix trimsignal = NULL;
static flappie_matrix normsignal = NULL;
static float * normsig_arr = NULL;



/**  Initialise test
 *
 *   @returns 0 on success, non-zero on failure
 **/
int init_test_signal(void) {
    rawsignal = read_flappie_matrix(rawsignalfile);
    if(NULL == rawsignal){
        return 1;
    }

    trimsignal = read_flappie_matrix(trimsignalfile);
    if(NULL == trimsignal){
        return 1;
    }

    normsignal = read_flappie_matrix(normsignalfile);
    if(NULL == normsignal){
        return 1;
    }
    normsig_arr = array_from_flappie_matrix(normsignal);

    return 0;
}

/**  Clean up after test
 *
 *   @returns 0 on success, non-zero on failure
 **/
int clean_test_signal(void) {
    free(normsig_arr);
    normsignal = free_flappie_matrix(normsignal);
    trimsignal = free_flappie_matrix(trimsignal);
    return 0;
}

void test_trim_signal(void) {
    const int winlen = 100;
    raw_table rt = {0};
    rt.raw  = array_from_flappie_matrix(rawsignal);
    CU_ASSERT_PTR_NOT_NULL_FATAL(rt.raw);
    rt.n = rt.end = rawsignal->nc;

    {   // Scale raw data to pA
        const float range = 1373.41f;
        const float digitisation = 8192.0f;
        const float unit = range / digitisation;
        const float offset = 16.0f;

        for(size_t i=0 ; i < rt.n ; i++){
            rt.raw[i] = (rt.raw[i] + offset) * unit;
        }
    }

    rt = trim_raw_by_mad(rt, winlen, 0.0f);
    CU_ASSERT_EQUAL(rt.start, 0);
    CU_ASSERT_EQUAL(rt.end, (rt.n / winlen) * winlen);

    rt.start += 200;
    rt.end -= 10;

    flappie_matrix mat_trim = mat_from_array(rt.raw + rt.start, 1, rt.end - rt.start);
    CU_ASSERT_PTR_NOT_NULL_FATAL(mat_trim);
    CU_ASSERT_EQUAL(mat_trim->nc, trimsignal->nc);

    CU_ASSERT_TRUE(equality_flappie_matrix(mat_trim, trimsignal, 1e-4));

    mat_trim = free_flappie_matrix(mat_trim);
    free(rt.raw);
}

void test_normalise_signal(void) {
    float * sigarr = array_from_flappie_matrix(trimsignal);
    size_t n = trimsignal->nc;
    CU_ASSERT_PTR_NOT_NULL_FATAL(sigarr);

    medmad_normalise_array(sigarr, n);
    CU_ASSERT_TRUE(equality_arrayf(sigarr, normsig_arr, n, 1e-5));

    free(sigarr);
}

static test_with_description tests[] = {
    {"Normalise trimmed signal", test_normalise_signal},
    {"Trimming of raw signal", test_trim_signal},
    {0}};

/**   Register tests with CUnit
 *
 *    @returns 0 on success, non-zero on failure
 **/
int register_test_signal(void) {
    return flappie_register_test_suite("Test manipulating signal", init_test_signal, clean_test_signal, tests);
}
