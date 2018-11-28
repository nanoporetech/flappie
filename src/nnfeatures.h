/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. License, v. 1.0. If a copy of the License was not
 *  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#pragma once
#ifndef FEATURES_H
#    define FEATURES_H
#    include <stdbool.h>
#    include "flappie_structures.h"
#    include "flappie_matrix.h"

flappie_matrix features_from_raw(const raw_table signal);

#endif /* FEATURES_H */
