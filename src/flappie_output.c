/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#include <err.h>
#include <stdlib.h>
#include <string.h>

#include "flappie_output.h"


enum flappie_outformat_type get_outformat(const char * formatstr){
    if(0 == strcmp(formatstr, "fasta")){
        return FLAPPIE_OUTFORMAT_FASTA;
    }
    if(0 == strcmp(formatstr, "fastq")){
        return FLAPPIE_OUTFORMAT_FASTQ;
    }
    if(0 == strcmp(formatstr, "sam")){
        return FLAPPIE_OUTFORMAT_SAM;
    }
    return FLAPPIE_OUTFORMAT_INVALID;
}


const char * flappie_outformat_string(enum flappie_outformat_type format){
    switch(format){
    case FLAPPIE_OUTFORMAT_FASTA:
        return "fasta";
    case FLAPPIE_OUTFORMAT_FASTQ:
        return "fastq";
    case FLAPPIE_OUTFORMAT_SAM:
        return "sam";
    case FLAPPIE_OUTFORMAT_INVALID:
        errx(EXIT_FAILURE, "Invalid flappie output %s:%d", __FILE__, __LINE__);
    default:
        errx(EXIT_FAILURE, "Flappie enum failure -- report bug\n");
    }

   return NULL;
}


void fprintf_format(enum flappie_outformat_type outformat, FILE * fp,
                    const char * uuid, const char *readname,
                    bool uuid_primary, const char * prefix,
                    const struct _raw_basecall_info res){
    switch(outformat){
    case FLAPPIE_OUTFORMAT_FASTA:
        fprintf_fasta(fp, uuid, readname, uuid_primary, prefix, res);
        break;
    case FLAPPIE_OUTFORMAT_FASTQ:
        fprintf_fastq(fp, uuid, readname, uuid_primary, prefix, res);
        break;
    case FLAPPIE_OUTFORMAT_SAM:
        fprintf_sam(fp, uuid, readname, uuid_primary, prefix, res);
        break;
    case FLAPPIE_OUTFORMAT_INVALID:
        errx(EXIT_FAILURE, "Invalid flappie output %s:%d", __FILE__, __LINE__);
    default:
        errx(EXIT_FAILURE, "Flappie enum failure -- report bug\n");
    }
}


void printf_format(enum flappie_outformat_type outformat,
		   const char * uuid, const char *readname,
		   bool uuid_primary, const char * prefix,
		   const struct _raw_basecall_info res){
    fprintf_format(outformat, stdout, uuid, readname, uuid_primary, prefix, res);
}


void fprint_string(FILE * fp, const char * str, bool newline){
    if(NULL ==fp || NULL == str){
        return;
    }
    while('\0' != *str){
        fputc(*str, fp);
        str++;
    }
    if(newline){
        fputc('\n', fp);
    }
}


void fprintf_fasta(FILE * fp, const char * uuid, const char *readname,
                   bool uuid_primary, const char * prefix,
                   const struct _raw_basecall_info res) {
    fprintf(fp, ">%s%s  { \"filename\" : \"%s\", \"uuid\" : \"%s\", \"normalised_score\" : %f,  \"nblock\" : %zu,  \"sequence_length\" : %zu,  \"blocks_per_base\" : %f, \"nsample\" : %zu, \"trim\" : [ %zu, %zu ] }\n",
                   prefix, uuid_primary ? uuid : readname, readname, uuid, 
                   -res.score / res.nblock, res.nblock, res.basecall_length,
                   (float)res.nblock / (float)res.basecall_length,
                   res.rt.n, res.rt.start, res.rt.end);
    fprint_string(fp, res.basecall, true);
    fflush(fp);
}


void fprintf_fastq(FILE * fp, const char * uuid, const char *readname,
                   bool uuid_primary, const char * prefix,
                   const struct _raw_basecall_info res) {
    if(NULL == res.quality){
        warnx("Can't output fastq for reads without quality values");
        return;
    }
    fprintf(fp, "@%s%s  { \"filename\" : \"%s\", \"uuid\" : \"%s\", \"normalised_score\" : %f,  \"nblock\" : %zu,  \"sequence_length\" : %zu,  \"blocks_per_base\" : %f, \"nsample\" : %zu, \"trim\" : [ %zu, %zu ] }\n",
                   prefix, uuid_primary ? uuid : readname, readname, uuid, 
                   -res.score / res.nblock, res.nblock, res.basecall_length,
                   (float)res.nblock / (float)res.basecall_length,
                   res.rt.n, res.rt.start, res.rt.end);
    fprint_string(fp, res.basecall, true);
    fputs("+\n", fp);
    fprint_string(fp, res.quality, true);
    fflush(fp);
}


void fprintf_sam(FILE * fp,  const char * uuid, const char *readname,
                 bool uuid_primary, const char * prefix,
                 const struct _raw_basecall_info res) {
    fprintf(fp, "%s%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", prefix,
                   uuid_primary ? uuid : readname, res.basecall, res.quality ? res.quality : "");
    fprint_string(fp, res.basecall, false);
    fputc('\t', fp);
    fprint_string(fp, res.quality, true);
    fflush(fp);
}
