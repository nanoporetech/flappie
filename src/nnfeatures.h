#pragma once
#ifndef FEATURES_H
#    define FEATURES_H
#    include <stdbool.h>
#    include "flappie_structures.h"
#    include "flappie_matrix.h"

flappie_matrix features_from_raw(const raw_table signal);

#endif /* FEATURES_H */
