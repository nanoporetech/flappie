/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#define BANANA 1
#include <CUnit/Basic.h>
#include <stdbool.h>

#include <util.h>
#include <test_common.h>

/**  Initialise test
 *
 *   @returns 0 on success, non-zero on failure
 **/
int init_test_util(void) {
    return 0;
}

/**  Clean up after test
 *
 *   @returns 0 on success, non-zero on failure
 **/
int clean_test_util(void) {
    return 0;
}

void test_median_odd_util(void) {
    float arr[5] = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f};
    float med = medianf(arr, 5);
    CU_ASSERT_EQUAL(med, 2.f);
}

void test_median_even_util(void) {
    float arr[5] = {0.0f, 1.0f, 2.0f, 3.0f};
    float med = medianf(arr, 4);
    CU_ASSERT_DOUBLE_EQUAL(med, 1.5f, 1e-5);
}

static test_with_description tests[] = {
    {"Median of odd length array", test_median_odd_util},
    {"Median of even length array", test_median_even_util},
    {0}};

/**   Register tests with CUnit
 *
 *    @returns 0 on success, non-zero on failure
 **/
int register_test_util(void) {
    return flappie_register_test_suite("Miscellaneous utility functions", init_test_util, clean_test_util, tests);
    return 0;
}
