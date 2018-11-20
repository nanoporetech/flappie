#include <dirent.h>
#include <fcntl.h>
#include <glob.h>
#include <libgen.h>
#include <math.h>
#include <stdio.h>
#include <strings.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "decode.h"
#include "fast5_interface.h"
#include "layers.h"
#include "networks.h"
#include "flappie_common.h"
#include "flappie_licence.h"
#include "flappie_output.h"
#include "flappie_stdlib.h"
#include "util.h"
#include "version.h"

#if !defined(FLAPPIE_VERSION)
#    define FLAPPIE_VERSION "unknown"
#endif
const char *argp_program_version = "flappie " FLAPPIE_VERSION;
const char *argp_program_bug_address = "<tim.massingham@nanoporetech.com>";

// Doesn't play nice with other headers, include last
#include <argp.h>


extern const char *argp_program_version;
extern const char *argp_program_bug_address;
static char doc[] = "Flappie basecaller -- basecall from raw signal";
static char args_doc[] = "fast5 [fast5 ...]";
static struct argp_option options[] = {
    //{"dump", 4, "filename", 0, "Dump annotated blocks to HDF5 file"},
    {"format", 'f', "format", 0, "Format to output reads (FASTA or SAM)"},
    {"limit", 'l', "nreads", 0, "Maximum number of reads to call (0 is unlimited)"},
    {"model", 'm', "name", 0, "Model to use (\"help\" to list)"},
    {"output", 'o', "filename", 0, "Write to file rather than stdout"},
    {"prefix", 'p', "string", 0, "Prefix to append to name of each read"},
    {"temperature", 7, "factor", 0, "Temperature for weights"},
    {"trim", 't', "start:end", 0, "Number of samples to trim, as start:end"},
    {"licence", 10, 0, 0, "Print licensing information"},
    {"license", 11, 0, OPTION_ALIAS, "Print licensing information"},
    {"segmentation", 3, "chunk:percentile", 0, "Chunk size and percentile for variance based segmentation"},
    {"hdf5-compression", 12, "level", 0,
     "Gzip compression level for HDF5 output (0:off, 1: quickest, 9: best)"},
    {"hdf5-chunk", 13, "size", 0, "Chunk size for HDF5 output"},

    {"uuid", 14, 0, 0, "Output UUID"},
    {"no-uuid", 15, 0, OPTION_ALIAS, "Output read file"},
    {0}
};



struct arguments {
    int compression_level;
    int compression_chunk_size;
    char * dump;
    enum flappie_outformat_type outformat;
    int limit;
    enum model_type model;
    FILE * output;
    char * prefix;
    float temperature;
    int trim_start;
    int trim_end;
    int varseg_chunk;
    float varseg_thresh;
    char ** files;
    bool uuid;
};

static struct arguments args = {
    .compression_level = 1,
    .compression_chunk_size = 200,
    .dump = NULL,
    .limit = 0,
    .model = FLAPPIE_MODEL_R941_NATIVE,
    .output = NULL,
    .outformat = FLAPPIE_OUTFORMAT_FASTQ,
    .prefix = "",
    .temperature = 1.0f,
    .trim_start = 200,
    .trim_end = 10,
    .varseg_chunk = 100,
    .varseg_thresh = 0.0f,
    .files = NULL,
    .uuid = false
};

void fprint_flappie_models(FILE * fh){
    if(NULL == fh){
        return;
    }

    for(size_t mdl=0 ; mdl < flappie_nmodel ; mdl++){
        fprintf(fh, "%10s : %s\n", flappie_model_string(mdl), flappie_model_description(mdl));
    }
}

static error_t parse_arg(int key, char * arg, struct  argp_state * state){
    int ret = 0;
    char * next_tok = NULL;

    switch(key){
    case 'd':
        args.dump = arg;
        break;
    case 'f':
        args.outformat = get_outformat(arg);
        if(FLAPPIE_OUTFORMAT_INVALID == args.outformat){
            errx(EXIT_FAILURE, "Unrecognised output format \"%s\".", arg);
        }
        break;
    case 'l':
        args.limit = atoi(arg);
        assert(args.limit > 0);
        break;
    case 'm':
        if(0 == strcasecmp(arg, "help")){
            fprint_flappie_models(stdout);
            exit(EXIT_SUCCESS);
        }
        args.model = get_flappie_model_type(arg);
        if(FLAPPIE_MODEL_INVALID == args.model){
            fprintf(stdout, "Invalid Flappie model \"%s\".\n", arg);
            fprint_flappie_models(stdout);
            exit(EXIT_FAILURE);
        }
        break;
    case 'o':
        args.output = fopen(arg, "w");
        if(NULL == args.output){
            errx(EXIT_FAILURE, "Failed to open \"%s\" for output.", arg);
        }
        break;
    case 'p':
        args.prefix = arg;
        break;
    case 't':
        args.trim_start = atoi(strtok(arg, ":"));
        next_tok = strtok(NULL, ":");
        if(NULL != next_tok){
            args.trim_end = atoi(next_tok);
        } else {
            args.trim_end = args.trim_start;
        }
        assert(args.trim_start >= 0);
        assert(args.trim_end >= 0);
        break;
    case 3:
        args.varseg_chunk = atoi(strtok(arg, ":"));
        next_tok = strtok(NULL, ":");
        if(NULL == next_tok){
            errx(EXIT_FAILURE, "--segmentation should be of form chunk:percentile");
        }
        args.varseg_thresh = atof(next_tok) / 100.0;
        assert(args.varseg_chunk >= 0);
        assert(args.varseg_thresh > 0.0 && args.varseg_thresh < 1.0);
        break;
    case 7:
	args.temperature = atof(arg);
	assert(isfinite(args.temperature) && args.temperature > 0.0f);
        break;
    case 10:
    case 11:
        ret = fputs(flappie_licence_text, stdout);
        exit((EOF != ret) ? EXIT_SUCCESS : EXIT_FAILURE);
        break;
    case 12:
        args.compression_level = atoi(arg);
        assert(args.compression_level >= 0 && args.compression_level <= 9);
        break;
    case 13:
        args.compression_chunk_size = atoi(arg);
        assert(args.compression_chunk_size > 0);
        break;
    case 14:
        args.uuid = true;
        break;
    case 15:
        args.uuid = false;
        break;
    case ARGP_KEY_NO_ARGS:
        argp_usage (state);
        break;

    case ARGP_KEY_ARG:
        args.files = &state->argv[state->next - 1];
        state->next = state->argc;
        break;

    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}


static struct argp argp = {options, parse_arg, args_doc, doc};


static struct _raw_basecall_info calculate_post(char * filename, enum model_type model){
    RETURN_NULL_IF(NULL == filename, (struct _raw_basecall_info){0});

    raw_table rt = read_raw(filename, true);
    RETURN_NULL_IF(NULL == rt.raw, (struct _raw_basecall_info){0});

    rt = trim_and_segment_raw(rt, args.trim_start, args.trim_end, args.varseg_chunk, args.varseg_thresh);
    RETURN_NULL_IF(NULL == rt.raw, (struct _raw_basecall_info){0});

    medmad_normalise_array(rt.raw + rt.start, rt.end - rt.start);

    flappie_matrix trans_weights = flipflop_transitions(rt, args.temperature, model);
    if (NULL == trans_weights) {
        free(rt.raw);
        free(rt.uuid);
        return (struct _raw_basecall_info){0};
    }

    const size_t nbase = nbase_from_flipflop_nparam(trans_weights->nr);
    const int nblock = trans_weights->nc;
    int * path = calloc(nblock + 2, sizeof(int));
    int * path_idx = calloc(nblock + 2, sizeof(int));
    float * qpath = calloc(nblock + 2, sizeof(float));
    int * pos = calloc(nblock + 1, sizeof(int));

    float score = NAN;

    flappie_matrix posterior = transpost_crf_flipflop(trans_weights, true);
    score = decode_crf_flipflop(posterior, false, path, qpath);
    size_t path_nidx = change_positions(path, nblock, path_idx);

    char * basecall = calloc(path_nidx + 1, sizeof(char));
    char * quality = calloc(path_nidx + 1, sizeof(char));
    for(size_t i=0 ; i < path_nidx ; i++){
        const size_t idx = path_idx[i];
        basecall[i] = base_lookup[path[idx] % nbase];
        quality[i] = phredf(expf(qpath[idx]));
    }

    posterior = free_flappie_matrix(posterior);
    free(qpath);
    free(path_idx);
    free(path);
    trans_weights = free_flappie_matrix(trans_weights);
    const size_t basecall_length = strlen(basecall);

    return (struct _raw_basecall_info) {
    	.score = score, 
        .rt = rt,
        .basecall = basecall,
        .quality = quality,
        .basecall_length = basecall_length,
        .pos = pos,
        .nblock = nblock};
}


int main(int argc, char * argv[]){
    argp_parse(&argp, argc, argv, 0, 0, NULL);
    if(NULL == args.output){
        args.output = stdout;
    }

    hid_t hdf5out = -1;
    if (NULL != args.dump) {
        int fd = open(args.dump, O_CREAT | O_WRONLY | O_EXCL, S_IRUSR | S_IWUSR);
        if(fd < 0){
            hdf5out = H5Fopen(args.dump, H5F_ACC_RDWR, H5P_DEFAULT);
        } else {
            close(fd);
            unlink(args.dump);
            hdf5out = H5Fcreate(args.dump, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        }
        if (hdf5out < 0){
            warnx("Failed to create \"%s\" for dumping.\n", args.dump);
            args.dump = NULL;
        }
    }


    int nfile = 0;
    for( ; args.files[nfile] ; nfile++);

    int reads_started = 0;
    const int reads_limit = args.limit;

    for(int fn=0 ; fn < nfile ; fn++){
        if(reads_limit > 0 && reads_started >= reads_limit){
            continue;
        }
        //  Iterate through all files and directories on command line.
        glob_t globbuf;
        {
            // Find all files matching commandline argument using system glob
            const size_t rootlen = strlen(args.files[fn]);
            char * globpath = calloc(rootlen + 9, sizeof(char));
            memcpy(globpath, args.files[fn], rootlen * sizeof(char));
            {
                DIR * dirp = opendir(args.files[fn]);
                if(NULL != dirp){
                    // If filename is a directory, add wildcard to find all fast5 files within it
                    memcpy(globpath + rootlen, "/*.fast5", 8 * sizeof(char));
                    closedir(dirp);
                }
            }
            int globret = glob(globpath, GLOB_NOSORT, NULL, &globbuf);
            free(globpath);
            if(0 != globret){
                if(GLOB_NOMATCH == globret){
                    warnx("File or directory \"%s\" does not exist or no fast5 files found.", args.files[fn]);
                }
                globfree(&globbuf);
                continue;
            }
        }

        for(size_t fn2=0 ; fn2 < globbuf.gl_pathc ; fn2++){
            if(reads_limit > 0 && reads_started >= reads_limit){
                continue;
            }
            reads_started += 1;

            char * filename = globbuf.gl_pathv[fn2];
            struct _raw_basecall_info res = calculate_post(filename, args.model);
            if(NULL == res.basecall){
                warnx("No basecall returned for %s", filename);
                continue;
            }

            fprintf_format(args.outformat, args.output, res.rt.uuid, 
                           basename(filename), args.uuid, args.prefix, res);

            write_trace(hdf5out, basename(filename), res,
                        args.compression_chunk_size,
                        args.compression_level);


            free(res.rt.raw);
            free(res.rt.uuid);
            free(res.basecall);
            free(res.quality);
            free(res.pos);
        }
        globfree(&globbuf);
    }


    if (hdf5out >= 0) {
        H5Fclose(hdf5out);
    }

    if(stdout != args.output){
        fclose(args.output);
    }

    return EXIT_SUCCESS;
}
