#pragma once
#ifndef FAST5_INTERFACE_H
#define FAST5_INTERFACE_H

#include <hdf5.h>
#include <stdbool.h>

#include "flappie_structures.h"

raw_table read_raw(const char *filename, bool scale_to_pA);
hid_t open_or_create_hdf5(const char * filename);

void write_summary(hid_t hdf5file, const char *readname,
                   const struct _raw_basecall_info res, 
                   hsize_t chunk_size, int compression_level);


#endif /* FAST5_INTERFACE_H */
