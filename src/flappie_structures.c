/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

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
