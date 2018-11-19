#pragma once
#ifndef TEST_COMMON_H
#    define TEST_COMMON_H

typedef struct {
    char * description;
    void (*testfunc)(void);
} test_with_description;

static int flappie_register_test_suite(char const *suite_name, int (*init_test)(void), int(*clean_test)(void), const test_with_description test[]){

    CU_pSuite suite = CU_add_suite(suite_name, init_test, clean_test);
    if (NULL == suite) {
        return CU_get_error();
    }

    if(NULL != test){
        for(size_t i=0 ; NULL != test[i].testfunc ; i++){
            if(NULL == CU_add_test(suite, test[i].description, test[i].testfunc)) {
                return CU_get_error();
            }
        }
    }

    return 0;
}
#endif /* TEST_COMMON_H */
