#pragma once
#ifndef NETWORKS_H
#    define NETWORKS_H
#    include <stdbool.h>
#    include "flappie_matrix.h"
#    include "flappie_structures.h"

flappie_matrix flipflop_transitions_r94pcr(const raw_table signal, float temperature);
flappie_matrix flipflop_transitions_r941native5mC(const raw_table signal, float temperature);


#endif    /* NETWORKS_H */
