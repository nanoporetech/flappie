#pragma once
#ifndef FAST5_INTERFACE_H
#define FAST5_INTERFACE_H

#include <hdf5.h>
#include <stdbool.h>

#include "flappie_structures.h"

raw_table read_raw(const char *filename, bool scale_to_pA);

void write_trace(hid_t hdf5file, const char *readname,
                 const struct _raw_basecall_info res, 
                 hsize_t chunk_size, int compression_level);


#endif /* FAST5_INTERFACE_H */
