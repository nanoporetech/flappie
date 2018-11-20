#include <stdlib.h>

#include "flappie_structures.h"

void free_raw_table(raw_table * tbl){
    free(tbl->uuid);
    free(tbl->raw);
}

void free_raw_basecall_info(struct _raw_basecall_info * ptr){
    free_raw_table(&(ptr->rt));
    free(ptr->basecall);
    free(ptr->quality);
    free(ptr->pos);
    free_flappie_imatrix(ptr->trace);
}
