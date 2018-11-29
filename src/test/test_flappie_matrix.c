/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#define BANANA 1
#include <CUnit/Basic.h>
#include <stdbool.h>

#include <flappie_matrix.h>
#include <test_common.h>

/**  Initialise test
 *
 *   @returns 0 on success, non-zero on failure
 **/
int init_test_flappie_matrix(void) {
    return 0;
}

/**  Clean up after test
 *
 *   @returns 0 on success, non-zero on failure
 **/
int clean_test_flappie_matrix(void) {
    return 0;
}

void test_rownormalise_flappie_matrix_helper(int nr) {
    flappie_matrix mat = make_flappie_matrix(nr, 1);
    CU_ASSERT_PTR_NOT_NULL_FATAL(mat);
    const int stride = mat->stride;
    for(int i=0 ; i < stride ; i++){
        mat->data.f[i] = 1.0f;
    }
    row_normalise_inplace(mat);

    const float expected = 1.0f / nr;
    for(int i=0 ; i < nr ; i++){
        CU_ASSERT_DOUBLE_EQUAL(mat->data.f[i], expected, 1e-5);
    }
    mat = free_flappie_matrix(mat);
}

void test_rownormalise_nr08flappie_matrix(void){
    test_rownormalise_flappie_matrix_helper(8);
}
void test_rownormalise_nr09flappie_matrix(void){
    test_rownormalise_flappie_matrix_helper(9);
}
void test_rownormalise_nr10flappie_matrix(void){
    test_rownormalise_flappie_matrix_helper(10);
}
void test_rownormalise_nr11flappie_matrix(void){
    test_rownormalise_flappie_matrix_helper(11);
}

static test_with_description tests[] = {
    {"Row normalisation edge case nr  8", test_rownormalise_nr08flappie_matrix},
    {"Row normalisation edge case nr  9", test_rownormalise_nr09flappie_matrix},
    {"Row normalisation edge case nr 10", test_rownormalise_nr10flappie_matrix},
    {"Row normalisation edge case nr 11", test_rownormalise_nr11flappie_matrix},
    {0}};

/**   Register tests with CUnit
 *
 *    @returns 0 on success, non-zero on failure
 **/
int register_test_matrix(void) {
    return flappie_register_test_suite("Functions on flappie matrices", init_test_flappie_matrix, clean_test_flappie_matrix, tests);
    return 0;
}
