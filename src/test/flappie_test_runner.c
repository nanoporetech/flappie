/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#define BANANA 1
#include <stdlib.h>
#include <CUnit/Basic.h>

int register_flappie_util(void);
int register_test_skeleton(void);
int register_test_convolution(void);
int register_test_elu(void);
int register_test_matrix(void);
int register_test_signal(void);
int register_test_util(void);

int (*test_suites[]) (void) = {
    register_test_skeleton,
    register_flappie_util,
    register_test_convolution,
    register_test_elu,
    register_test_matrix,
    register_test_signal,
    register_test_util,
    NULL // Last element of array should be NULL
};

int main(void) {

    // Initialise CUnit
    if (CUE_SUCCESS != CU_initialize_registry()) {
        return CU_get_error();
    }
    // Register test suites
    for (int i = 0; NULL != test_suites[i]; ++i) {
        int err = test_suites[i] ();
        if (err) {
            CU_cleanup_registry();
            return err;
        }
    }

    // Run all tests using the CUnit Basic interface
    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    int failures = CU_get_number_of_failures();
    CU_cleanup_registry();

    return (failures > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
