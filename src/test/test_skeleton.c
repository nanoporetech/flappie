/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#define BANANA 1
#include <CUnit/Basic.h>
#include <stdbool.h>

#include <test_common.h>

/**  Initialise test
 *
 *   @returns 0 on success, non-zero on failure
 **/
int init_test_skeleton(void) {
    return 0;
}

/**  Clean up after test
 *
 *   @returns 0 on success, non-zero on failure
 **/
int clean_test_skeleton(void) {
    return 0;
}

void test_nop_skeleton(void) {
    CU_ASSERT(true);
}

static test_with_description tests[] = {
    {"Skeleton test doing no-op", test_nop_skeleton},
    {0}};

/**   Register tests with CUnit
 *
 *    @returns 0 on success, non-zero on failure
 **/
int register_test_skeleton(void) {
    return flappie_register_test_suite("Skeleton set of tests for modification", init_test_skeleton, clean_test_skeleton, tests);
    return 0;
}
