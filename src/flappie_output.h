/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#pragma once
#ifndef FLAPPIE_OUTPUT_H
#define FLAPPIE_OUTPUT_H

#include <stdbool.h>
#include <stdio.h>

#include "flappie_structures.h"

enum flappie_outformat_type {FLAPPIE_OUTFORMAT_FASTA,
                             FLAPPIE_OUTFORMAT_FASTQ,
                             FLAPPIE_OUTFORMAT_SAM,
                             FLAPPIE_OUTFORMAT_INVALID};

enum flappie_outformat_type get_outformat(const char * formatstr);
const char * flappie_outformat_string(enum flappie_outformat_type format);


void printf_format(enum flappie_outformat_type outformat,
                   const char * uuid, const char *readname,
                   bool uuid_primary, const char * prefix,
                   const struct _raw_basecall_info res);

void fprintf_format(enum flappie_outformat_type outformat, FILE * fp,
                    const char * uuid, const char *readname,
                    bool uuid_primary, const char * prefix,
                    const struct _raw_basecall_info res);

void fprintf_fasta(FILE * fp, const char * uuid, const char *readname,
                   bool uuid_primary, const char * prefix,
                   const struct _raw_basecall_info res);

void fprintf_fastq(FILE * fp, const char * uuid, const char *readname,
                   bool uuid_primary, const char * prefix,
                   const struct _raw_basecall_info res);

void fprintf_sam(FILE * fp,  const char * uuid, const char *readname,
                 bool uuid_primary, const char * prefix,
                 const struct _raw_basecall_info res);

#endif /* FLAPPIE_OUTPUT_H */
