/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#include <stdio.h>

#include "decode.h"
#include "layers.h"
#include "flappie_stdlib.h"
#include "util.h"


float argmax_decoder(const_flappie_matrix logpost, int *seq) {
    RETURN_NULL_IF(NULL == logpost, NAN);
    RETURN_NULL_IF(NULL == seq, NAN);

    const size_t nblock = logpost->nc;
    const size_t nstate = logpost->nr;
    assert(nstate > 0);
    const size_t stride = logpost->stride;
    assert(stride > 0);

    float logscore = 0;
    for (size_t blk = 0; blk < nblock; blk++) {
        const size_t offset = blk * stride;
        const int imax = argmaxf(logpost->data.f + offset, nstate);
        logscore += logpost->data.f[offset + imax];
        seq[blk] = (imax == nstate - 1) ? -1 : imax;
    }

    return logscore;
}


char * collapse_repeats(int const * path, size_t npos, int modbase){
    assert(modbase > 0);
    RETURN_NULL_IF(NULL == path, NULL);

    int nbase = 1;
    for(size_t pos=1 ; pos < npos ; pos++){
        if(path[pos] != path[pos - 1]){
            nbase += 1;
        }
    }

    char * basecall = calloc(nbase + 1, sizeof(char));
    RETURN_NULL_IF(NULL == basecall, NULL);

    basecall[0] = base_lookup[path[0] % modbase];
    for(size_t pos=1, bpos=1 ; pos < npos ; pos++){
        if(path[pos] != path[pos - 1]){
            assert(bpos < nbase);
            basecall[bpos] = basechar(path[pos] % modbase);
            bpos += 1;
        }
    }

    return basecall;
}


size_t change_positions(int const * path, size_t npos, int * chpos){
    RETURN_NULL_IF(NULL == path, 0);
    RETURN_NULL_IF(NULL == chpos, 0);

    size_t nch = 0;
    for(size_t pos=1 ; pos < npos ; pos++){
        if(path[pos] == path[pos -1 ]){
            continue;
        }
        chpos[nch] = pos;
        nch += 1;
    }
    return nch;
}


void colmaxf(float * x, int nr, int nc, int * idx){
    assert(nr > 0);
    assert(nc > 0);
    RETURN_NULL_IF(NULL == x,);
    RETURN_NULL_IF(NULL == idx,);

    for(int r=0 ; r < nr ; r++){
        // Initialise
        idx[r] = 0;
    }

    for(int c=1 ; c < nc ; c++){
        const size_t offset2 = c * nr;
        for(int r=0 ; r<nr ; r++){
            if(x[offset2 + r] > x[idx[r] * nr + r]){
                idx[r] = c;
            }
        }
    }
}


inline size_t trans_lookup(size_t from, size_t to, size_t nbase){
    assert(nbase >= 0);
    assert(from >= 0 && from < nbase + nbase);
    assert(to >= 0 && to < nbase + nbase);
    assert(to < nbase || ((to % nbase) == (from % nbase)));

    const size_t nstate = nbase + nbase;
    const size_t offset = nbase * nstate;

    return (to < nbase) ? (to * nstate + from) : (offset + from);
}


/**   Viterbi decoding of CRF flipflop
 **/
float decode_crf_flipflop(const_flappie_matrix trans, bool combine_stays, int * path, float * qpath){
    RETURN_NULL_IF(NULL == trans, NAN);
    RETURN_NULL_IF(NULL == path, NAN);
    RETURN_NULL_IF(NULL == qpath, NAN);

    const size_t nblk = trans->nc;
    const size_t nbase = roundf((-1.0f + sqrtf(1.0f + 2.0f * trans->nr)) / 2.0f);
    const size_t nstate = nbase + nbase;
    assert(nstate == nbase + nbase);
    assert(nstate * (nbase + 1) == trans->nr);

    float * mem = calloc(2 * nstate, sizeof(float));
    flappie_imatrix tb = make_flappie_imatrix(nstate, nblk);
    if(NULL == mem || NULL == tb){
        tb = free_flappie_imatrix(tb);
        free(mem);
        return NAN;
    }

    float * curr = mem;
    float * prev = mem + nstate;


    //  Forwards Viterbi pass
    for(size_t blk=0 ; blk < nblk ; blk++){
        const size_t offset = blk * trans->stride;
        const size_t offset_flop = offset + nstate * nbase;
        const size_t tboffset = blk * tb->stride;
        {   // Swap
            float * tmp = curr;
            curr = prev;
            prev = tmp;
        }

        for(size_t b2=nbase ; b2 < nstate ; b2++){
            // Stay in flop state
            curr[b2] = prev[b2] + trans->data.f[offset_flop + b2];
            tb->data.f[tboffset + b2] = b2;
            // Move from flip to flop state
            const size_t from_base = b2 - nbase;
            const float score = prev[from_base] + trans->data.f[offset_flop + from_base];
            if(score > curr[b2]){
                curr[b2] = score;
                tb->data.f[tboffset + b2] = from_base;
            }
        }


        for(size_t b1=0 ; b1 < nbase ; b1++){
	    //   b1 -- flip state
            const size_t offset_state = offset + b1 * nstate;
            curr[b1] = trans->data.f[offset_state + 0] + prev[0];
            tb->data.f[tboffset + b1] = 0;
            for(size_t from_state=1 ; from_state < nstate ; from_state++){
                // from_state either flip or flop
                const float score = trans->data.f[offset_state + from_state] + prev[from_state];
                if(score > curr[b1]){
                    curr[b1] = score;
                    tb->data.f[tboffset + b1] = from_state;
                }
            }
        }
    }

    //  Traceback
    const float score = valmaxf(curr, nstate);
    path[nblk] = argmaxf(curr, nstate);
    for(size_t blk=nblk ; blk > 0 ; blk--){
        const size_t offset = (blk - 1) * tb->stride;
        const size_t qoffset = (blk - 1) * trans->stride;
        path[blk - 1] = tb->data.f[offset + path[blk]];
        qpath[blk] = trans->data.f[qoffset + trans_lookup(path[blk-1], path[blk], nbase)];
    }
    qpath[0] = NAN;

    if(combine_stays){
        for(size_t blk=0 ; blk <= nblk ; blk++){
            path[blk] = (path[blk] < nbase) ? path[blk] : -1;
        }
    }

    tb = free_flappie_imatrix(tb);
    free(mem);

    return score;
}


/**   Decoding of CRF flip-posteriors with transition constraint
 **/
float constrained_crf_flipflop(const_flappie_matrix post, int * path){
    RETURN_NULL_IF(NULL == post, NAN);
    RETURN_NULL_IF(NULL == path, NAN);

    const size_t nblk = post->nc;
    const size_t nstate = post->nr;
    const size_t nbase = nstate / 2;
    assert(nstate == nbase + nbase);

    float * mem = calloc(2 * nstate, sizeof(float));
    flappie_imatrix tb = make_flappie_imatrix(nstate, nblk);
    if(NULL == mem || NULL == tb){
        tb = free_flappie_imatrix(tb);
        free(mem);
        return NAN;
    }

    float * curr = mem;
    float * prev = mem + nstate;


    //  Forwards Viterbi pass
    for(size_t blk=0 ; blk < nblk ; blk++){
        const size_t offset = blk * post->stride;
        const size_t tboffset = blk * tb->stride;
        {   // Swap
            float * tmp = curr;
            curr = prev;
            prev = tmp;
        }

        for(size_t b2=nbase ; b2 < nstate ; b2++){
            const size_t best = (prev[b2] > prev[b2 - nbase]) ? b2 : (b2 - nbase);
            curr[b2] = prev[best];
            tb->data.f[tboffset + b2] = best;
        }


        for(size_t b1=0 ; b1 < nbase ; b1++){
            const int from_best = argmaxf(prev, nstate);
            curr[b1] = prev[from_best];
            tb->data.f[tboffset + b1] = from_best;
        }

        for(size_t st=0 ; st < nstate ; st++){
            curr[st] += post->data.f[offset + st];
        }
    }

    //  Traceback
    const float score = valmaxf(curr, nstate);
    path[nblk] = argmaxf(curr, nstate);
    for(size_t blk=nblk ; blk > 0 ; blk--){
        const size_t offset = (blk - 1) * tb->stride;
        path[blk - 1] = tb->data.f[offset + path[blk]];
    }

    tb = free_flappie_imatrix(tb);
    free(mem);

    return score;
}


/**   Posterior probabilities of CRF flipflop
 **/
flappie_matrix posterior_crf_flipflop(const_flappie_matrix trans, bool return_log){
    RETURN_NULL_IF(NULL == trans, NULL);

    const size_t nblk = trans->nc;
    const size_t nbase = roundf((-1.0f + sqrtf(1.0f + 2.0f * trans->nr)) / 2.0f);
    const size_t nstate = nbase + nbase;
    assert(nstate == nbase + nbase);
    assert(nstate * (nbase + 1) == trans->nr);

    flappie_matrix fwd = make_flappie_matrix(nstate, nblk + 1);
    RETURN_NULL_IF(NULL == fwd, NULL);


    //  Forwards pass
    for(size_t blk=0 ; blk < nblk ; blk++){
        const size_t offset = blk * trans->stride;
        const size_t offset_flop = offset + nstate * nbase;

        float * prev = fwd->data.f + blk * fwd->stride;
        float * curr = prev + fwd->stride;

        for(size_t b2=nbase ; b2 < nstate ; b2++){
            // Stay in flop state
            curr[b2] = prev[b2] + trans->data.f[offset_flop + b2];
            // Move from flip to flop state
            const size_t from_base = b2 - nbase;
            const float score = prev[from_base] + trans->data.f[offset_flop + from_base];
            curr[b2] = logsumexpf(curr[b2], score);
        }


        for(size_t b1=0 ; b1 < nbase ; b1++){
	    //   b1 -- flip state
            const size_t offset_state = offset + b1 * nstate;
            curr[b1] = trans->data.f[offset_state + 0] + prev[0];
            for(size_t from_state=1 ; from_state < nstate ; from_state++){
                // from_state either flip or flop
                const float score = trans->data.f[offset_state + from_state] + prev[from_state];
                curr[b1] = logsumexpf(curr[b1], score);
            }
        }
    }

    float * mem = calloc(2 * nstate, sizeof(float));
    if(NULL == mem){
        free(fwd);
        return NULL;
    }
    float * prev = mem;
    float * curr = mem + nstate;

    //  Backwards pass
    for(size_t blk=nblk ; blk > 0 ; blk--){
        const size_t foffset = (blk - 1) * fwd->stride;
        const size_t offset = (blk - 1) * trans->stride;
        const size_t offset_flop = offset + nstate * nbase;

        {  // Swap
           float * tmp = prev;
           prev = curr;
           curr = tmp;
        }

        for(size_t b2=nbase ; b2 < nstate ; b2++){
            const size_t from_base = b2 - nbase;
            // Stay in flop state
            curr[b2] = prev[b2] + trans->data.f[offset_flop + b2];
            // Move from flip to flop state
            curr[from_base] = prev[b2] + trans->data.f[offset_flop + from_base];
        }


        for(size_t b1=0 ; b1 < nbase ; b1++){
	    //   b1 -- flip state
            const size_t offset_state = offset + b1 * nstate;
            for(size_t from_state=0 ; from_state < nstate ; from_state++){
                // from_state either flip or flop
                const float score = trans->data.f[offset_state + from_state] + prev[b1];
                curr[from_state] = logsumexpf(curr[from_state], score);
            }
        }

        for(size_t st=0 ; st < nstate ; st++){
            // Add to fwd vector
            fwd->data.f[foffset + st] += curr[st];
        }
    }

    free(mem);

    if(!return_log){
        exp_activation_inplace(fwd);
        row_normalise_inplace(fwd);
    }


    return fwd;
}


/**   Posterior probabilities of CRF flipflop
 **/
flappie_matrix transpost_crf_flipflop(const_flappie_matrix trans, bool return_log){
    RETURN_NULL_IF(NULL == trans, NULL);

    const size_t nblk = trans->nc;
    const size_t nbase = roundf((-1.0f + sqrtf(1.0f + 2.0f * trans->nr)) / 2.0f);
    const size_t nstate = nbase + nbase;
    assert(nstate == nbase + nbase);
    assert(nstate * (nbase + 1) == trans->nr);

    flappie_matrix fwd = make_flappie_matrix(nstate, nblk + 1);
    flappie_matrix tpost = make_flappie_matrix(trans->nr, nblk);
    if(NULL == fwd || NULL == tpost){
        fwd = free_flappie_matrix(fwd);
        tpost = free_flappie_matrix(tpost);
        return NULL;
    }


    //  Forwards pass
    for(size_t blk=0 ; blk < nblk ; blk++){
        const size_t offset = blk * trans->stride;
        const size_t offset_flop = offset + nstate * nbase;

        float * prev = fwd->data.f + blk * fwd->stride;
        float * curr = prev + fwd->stride;

        for(size_t b2=nbase ; b2 < nstate ; b2++){
            // Stay in flop state
            curr[b2] = prev[b2] + trans->data.f[offset_flop + b2];
            // Move from flip to flop state
            const size_t from_base = b2 - nbase;
            const float score = prev[from_base] + trans->data.f[offset_flop + from_base];
            curr[b2] = logsumexpf(curr[b2], score);
        }


        for(size_t b1=0 ; b1 < nbase ; b1++){
	    //   b1 -- flip state
            const size_t offset_state = offset + b1 * nstate;
            curr[b1] = trans->data.f[offset_state + 0] + prev[0];
            for(size_t from_state=1 ; from_state < nstate ; from_state++){
                // from_state either flip or flop
                const float score = trans->data.f[offset_state + from_state] + prev[from_state];
                curr[b1] = logsumexpf(curr[b1], score);
            }
        }
    }

    float * mem = calloc(2 * nstate, sizeof(float));
    if(NULL == mem){
        free(fwd);
        return NULL;
    }
    float * prev = mem;
    float * curr = mem + nstate;

    //  Backwards pass
    for(size_t blk=nblk ; blk > 0 ; blk--){
        const size_t foffset = (blk - 1) * fwd->stride;
        const size_t offset = (blk - 1) * trans->stride;
        const size_t offset_flop = offset + nstate * nbase;

        {  // Swap
           float * tmp = prev;
           prev = curr;
           curr = tmp;
        }


        //  Create tpost
        for(size_t b1=0 ; b1 < nbase ; b1++){
            //  End up in flip state
            const size_t offset_state = offset + b1 * nstate;
            for(size_t st=0 ; st < nstate ; st++){
                tpost->data.f[offset_state + st] = fwd->data.f[foffset + st] + prev[b1]
                                                 + trans->data.f[offset_state + st];
            }
        }
        for(size_t b=nbase ; b < nstate ; b++){
            //  End up in flop state
            const size_t fb = b - nbase;
            tpost->data.f[offset_flop + b] = fwd->data.f[foffset + b] + prev[b]
                                           + trans->data.f[offset_flop + b];
            tpost->data.f[offset_flop + fb] = fwd->data.f[foffset + fb] + prev[b]
                                           + trans->data.f[offset_flop + fb];
        }


        //  Update backwards vector
        for(size_t b2=nbase ; b2 < nstate ; b2++){
            const size_t from_base = b2 - nbase;
            // Stay in flop state
            curr[b2] = prev[b2] + trans->data.f[offset_flop + b2];
            // Move from flip to flop state
            curr[from_base] = prev[b2] + trans->data.f[offset_flop + from_base];
        }


        for(size_t b1=0 ; b1 < nbase ; b1++){
	    //   b1 -- flip state
            const size_t offset_state = offset + b1 * nstate;
            for(size_t from_state=0 ; from_state < nstate ; from_state++){
                // from_state either flip or flop
                const float score = trans->data.f[offset_state + from_state] + prev[b1];
                curr[from_state] = logsumexpf(curr[from_state], score);
            }
        }
    }


    free(mem);
    fwd = free_flappie_matrix(fwd);


    log_row_normalise_inplace(tpost);
    if(!return_log){
        exp_activation_inplace(tpost);
    }

    return tpost;
}

flappie_imatrix trace_from_posterior(const flappie_matrix tpost){
    RETURN_NULL_IF(NULL == tpost, NULL);
    const size_t nbase = nbase_from_flipflop_nparam(tpost->nr);
    const size_t nstate = nbase + nbase;
    assert((nbase + 1) * nstate == tpost->nr);


    flappie_imatrix trace = make_flappie_imatrix(nstate, tpost->nc + 1);
    RETURN_NULL_IF(NULL == trace, NULL);


    //  First Position
    for(size_t st_from=0 ; st_from < nstate ; st_from++){
        float sum = 0.0f;
        for(size_t st_to=0 ; st_to < nbase ; st_to++){
            sum += tpost->data.f[st_to * nstate + st_from];
        }
        sum += tpost->data.f[nbase * nstate + st_from];
        trace->data.f[st_from] = roundf(255.0f * sum);
    }

    //  Other positions
    for(size_t blk=0 ; blk < tpost->nc ; blk++){
        const size_t offset_trace = (blk + 1) * trace->stride;
        const size_t offset_post = blk * tpost->stride;
        for(size_t st_to=0 ; st_to < nbase ; st_to++){
            //  Transition to flip state
            const size_t offset2 = offset_post + st_to * nstate;
            float sum = tpost->data.f[offset2];
            for(size_t st_from=1 ; st_from < nstate ; st_from++){
                sum += tpost->data.f[offset2 + st_from];
            }
            trace->data.f[offset_trace + st_to] = roundf(255.0f * sum);
        }

        const size_t offset_post2 = blk * tpost->stride + nbase * nstate;
        for(size_t st_to=nbase ; st_to < nstate ; st_to++){
            const float sum = tpost->data.f[offset_post2 + (st_to - nbase)]
                            + tpost->data.f[offset_post2 + st_to];;
            trace->data.f[offset_trace + st_to] = roundf(255.0f * sum);
        }
    }

    return trace;
}


/**  Approximate mean of a Discrete Weibull distribution
 *
 *   @param shape  Shape parameter of distribution
 *   @param scale  Scale parameter of distribution
 *   @param maxval Maximum value to calculate up to
 **/
float dwmean(float shape, float scale, int maxval){
    assert(shape > 0.0f);
    assert(scale > 0.0f);
    assert(maxval > 0);
    float m = 0.0f;
    for(int i=1 ; i <= maxval ; i++){
        m += expf(-powf((float)i / scale, shape));
    }
    return m;
}


/**  Calculate length of base runs given path
 *
 *   For each non-stay element on path, calculated expected length of run
 *
 *   @param param  Flappie matrix [13 x nblk] containing predicted parameters
 *   @param path Array[nblk] containing best path; -1 for stay
 *   @param runlength[out] Array[nblk] to write runs out to; 0 for stay
 *
 *   @returns score of best path
 **/
size_t runlengths_mean(const_flappie_matrix param, const int * path, int * runlength){
    RETURN_NULL_IF(NULL == param, 0);
    RETURN_NULL_IF(NULL == path, 0);
    RETURN_NULL_IF(NULL == runlength, 0);

    const size_t nblk = param->nc;
    const size_t nparam = param->nr;
    const size_t nbase = nbase_from_runlength_nparam(nparam);
    const size_t param_shape_offset = 0;
    const size_t param_scale_offset = nbase;

    memset(runlength, 0, nblk * sizeof(int));

    size_t seqlen = 0;
    for(size_t blk=0 ; blk < nblk ; blk++){
        if(path[blk] < 0){
            // Short circuit stays
            continue;
        }
        const size_t offset = blk * param->stride + path[blk];
        const float shape = param->data.f[offset + param_shape_offset];
        const float scale = param->data.f[offset + param_scale_offset];
        const float meanest = dwmean(shape, scale, 100);
        runlength[blk] = 1 + roundf(meanest);
        seqlen += runlength[blk];
    }
    return seqlen;
}


/**  Calculate length of base runs given path
 *
 *   For each non-stay element on path, assign unit run
 *
 *   @param param  Flappie matrix [13 x nblk] containing predicted parameters
 *   @param path Array[nblk] containing best path; -1 for stay
 *   @param runlength[out] Array[nblk] to write runs out to; 0 for stay
 *
 *   @returns score of best path
 **/
size_t runlengths_unit(const_flappie_matrix param, const int * path, int * runlength){
    RETURN_NULL_IF(NULL == param, 0);
    RETURN_NULL_IF(NULL == path, 0);
    RETURN_NULL_IF(NULL == runlength, 0);

    const size_t nblk = param->nc;
    memset(runlength, 0, nblk * sizeof(int));

    size_t seqlen = 0;
    for(size_t blk=0 ; blk < nblk ; blk++){
        if(path[blk] < 0){
            // Short circuit stays
            continue;
        }
        runlength[blk] = 1;
        seqlen += runlength[blk];
    }
    return seqlen;
}


/**  Convert path and runlength arrays into base-space sequence
 *
 *   @param path Array[nblk] containing path; -1 for stay
 *   @param runlength Array[nblk] containing runlengths; 0 for stay
 *   @param nblk Number of blocks (length of two arrays)
 *
 *   @returns Null terminated array containing base-space sequence
 **/
char * runlength_to_basecall(const int * path, const int * runlength, size_t nblk){
    RETURN_NULL_IF(NULL == path, NULL);
    RETURN_NULL_IF(NULL == runlength, NULL);


    int seqlen = 0;
    for(size_t blk=0 ; blk < nblk ; blk++){
        seqlen += runlength[blk];
    }

    char * seq = calloc(seqlen + 1, sizeof(char));
    RETURN_NULL_IF(NULL == seq, NULL);

    for(size_t blk=0, i=0 ; blk < nblk ; blk++){
        if(path[blk] < 0){
            continue;
        }
        const char base = base_lookup[path[blk]];
        for(size_t rl=0 ; rl < runlength[blk] ; rl++, i++){
            seq[i] = base;
        }
    }

    return seq;
}


/**  Decoding of runlength model with multiple stay states
 *
 *   Runlength models have the restriction that you can't move to the same
 *   base that you are already in. Once in a move state, the model must
 *   immediately exit into either the corresponding stay state or a new
 *   move state.
 *
 *   Since this model is globally normalised, sum(exp(weights)) may not sum
 *   to one.
 *
 *   Order of params:  Each column contained 16 entries in the following order:
 *       0 --  3 : 'R' parameter for negative binomial distribution [ACGT]
 *       4 --  7 : 'P' parameter for negative binomial distribution [ACGT]
 *       8 -- 11 : Move weights [ACGT]
 *      12 -- 15 : Stay weights [ACGT]
 *
 *   @param param  Flappie matrix [16 x nblk] containing predicted parameters
 *   @param combine_stays Combine all stay states in path to single state
 *   @param path[out] Array to write out best path
 *
 *   @returns score of best path
 **/
float decode_runlength(const_flappie_matrix param, int * path){
    RETURN_NULL_IF(NULL == param, NAN);
    RETURN_NULL_IF(NULL == path, NAN);
    const size_t nblk = param->nc;
    const size_t nparam = param->nr;
    const size_t nbase = nbase_from_runlength_nparam(nparam);

    float * mem = calloc(2 * nbase, sizeof(float));
    char * traceback = calloc(nbase * nblk, sizeof(char));
    if(NULL == mem || NULL == traceback){
        free(mem);
        free(traceback);
        return NAN;
    }

    float * prev = mem;
    float * curr = mem + nbase;

    for(size_t blk=0 ; blk < nblk ; blk++){
        const size_t offset_move = blk * param->stride + nbase + nbase;
        const size_t offset_stay = offset_move + nbase;
        const size_t toffset = blk * nbase;
        {   //  Swap memory
            float * tmp;
            tmp = prev;
            prev = curr;
            curr = tmp;
        }

        {   // Move to new base (need best and second best indexes: gross modification of prev)
            const int idx = argmaxf(prev, nbase);
            const float max_score = prev[idx];
            prev[idx] = -HUGE_VAL;
            const int idx2 = argmaxf(prev, nbase);
            prev[idx] = max_score;

            for(size_t b=0 ; b < nbase ; b++){
                curr[b] = max_score;
                traceback[toffset + b] = idx;
            }
            curr[idx] = prev[idx2];
            traceback[toffset + idx] = idx2;

            for(size_t b=0 ; b < nbase ; b++){
                //  Add on move score
                curr[b] += param->data.f[offset_move + b];
            }
        }
        for(int b=0 ; b < nbase ; b++){
            //  Stay state : come either from corresponding move state or same stay
            const float stay_score = prev[b] +  param->data.f[offset_stay + b];
            if(stay_score > curr[b]){
                // Come from stay
                curr[b] = stay_score;
                traceback[toffset + b] = b + nbase;
            }
        }
    }


    for(size_t blk=0 ; blk < nblk ; blk++){
        path[blk] = -1;
    }
    size_t last_state = argmaxf(curr, nbase);
    const float logscore = curr[last_state];
    for(size_t blk=nblk ; blk > 0 ; blk--){
        const size_t blkm1 = blk - 1;
        const char state = traceback[blkm1 * nbase + last_state];
        if(state < nbase){
            path[blkm1] = last_state;
            last_state = state;
        }
    }

    free(traceback);
    free(mem);

    return logscore;
}


/**  Posterior state probabilities runlength model with multiple stay states
 *
 *   Runlength models have the restriction that you can't move to the same
 *   base that you are already in. Once in a move state, the model must
 *   immediately exit into either the corresponding stay state or a new
 *   move state.
 *
 *   Since this model is globally normalised, sum(exp(weights)) may not sum
 *   to one.
 *
 *   Order of params:  Each column contained 16 entries in the following order:
 *       0 --  3 : 'R' parameter for negative binomial distribution [ACGT]
 *       4 --  7 : 'P' parameter for negative binomial distribution [ACGT]
 *       8 -- 11 : Move weights [ACGT]
 *      12 -- 15 : Stay weights [ACGT]
 *
 *   @param param  Flappie matrix [16 x nblk] containing predicted parameters
 *   @param combine_stays Combine all stay states in path to single state
 *   @param path[out] Array to write out best path
 *
 *   @returns Flappie matrix containing log posterior probabilities
 **/
flappie_matrix posterior_runlength(const_flappie_matrix param){
    RETURN_NULL_IF(NULL == param, NULL);
    const size_t nblk = param->nc;
    const size_t nparam = param->nr;
    const size_t nbase = nbase_from_runlength_nparam(nparam);
    const size_t param_p_offset = nbase;
    const size_t param_cat_offset = param_p_offset + nbase;
    const size_t param_stay_offset = param_cat_offset + nbase;

    flappie_matrix fwd = make_flappie_matrix(nbase, nblk + 1);
    flappie_matrix post = make_flappie_matrix(nparam, nblk + 1);
    float * mem = calloc(nbase + nbase, sizeof(float));
    if(NULL == fwd || NULL == post || NULL == mem){
        fwd = free_flappie_matrix(fwd);
        post = free_flappie_matrix(post);
        free(mem);
        return NULL;
    }


    for(size_t blk=0 ; blk < nblk ; blk++){
        //  Forwards calculation
        const size_t offset = blk * param->stride + param_cat_offset;
        const size_t offset_stay = blk * param->stride + param_stay_offset;
        float * prev = fwd->data.f + blk * fwd->stride;
        float * curr = prev + fwd->stride;

        for(size_t b1=0 ; b1 < nbase ; b1++){
            curr[b1] = -HUGE_VAL;
            // Move from different base
            for(size_t b2=0 ; b2 < nbase ; b2++){
                if(b1 != b2){
                    curr[b1] = logsumexpf(curr[b1], prev[b2]);
                }
            }
            curr[b1] += param->data.f[offset + b1];
        }
        for(size_t b=0 ; b < nbase ; b++){
            // Stay in same base
            curr[b] = logsumexpf(curr[b], prev[b] + param->data.f[offset_stay + b]);
        }
    }


    float * prev = mem;
    float * curr = prev + nbase;
    for(size_t blk=nblk ; blk > 0 ; blk--){
        const size_t foffset = (blk - 1) * fwd->stride;
        const size_t offset = (blk - 1) * param->stride + param_cat_offset;
        const size_t offset_stay = (blk - 1) * param->stride + param_stay_offset;

        // Backwards
        {  // swap
            float * tmp = curr;
            curr = prev;
            prev = tmp;
        }

        for(size_t b1=0 ; b1 < nbase ; b1++){
            curr[b1] = -HUGE_VAL;
            post->data.f[offset + b1] = -HUGE_VAL;
            // Move from different base
            for(size_t b2=0 ; b2 < nbase ; b2++){
                if(b1 != b2){
                    curr[b1] = logsumexpf(curr[b1], prev[b2] + param->data.f[offset + b2]);
                    post->data.f[offset + b1] = logsumexpf(post->data.f[offset + b1], fwd->data.f[foffset + b2]);
                }
            }
            post->data.f[offset + b1] += prev[b1] + param->data.f[offset + b1];
        }
        for(size_t b=0 ; b < nbase ; b++){
            // Stay in same base
            curr[b] = logsumexpf(curr[b], prev[b] + param->data.f[offset_stay + b]);
            post->data.f[offset_stay + b] = fwd->data.f[foffset + b] + param->data.f[offset_stay + b] + prev[b];
        }

        float score = curr[0] + fwd->data.f[foffset + 0];
        for(size_t i=1 ; i < nbase ; i++){
            score = logsumexpf(score, curr[i] + fwd->data.f[foffset + i]);
        }
        //printf("score for blk %zu : %f\n", blk - 1, score);
    }

    const size_t last_offset = nblk * fwd->stride;
    float scoreF = fwd->data.f[last_offset];
    float scoreB = curr[0];
    for(size_t st=1 ; st < nbase ; st ++){
        scoreF = logsumexpf(scoreF, fwd->data.f[last_offset + st]);
        scoreB = logsumexpf(scoreB, curr[st]);
    }
    //printf("ScoreF = %f  scoreB = %f\n", scoreF, scoreB);

    free(mem);
    fwd = free_flappie_matrix(fwd);


    return post;
}


static inline size_t rle_trans_lookup(size_t base_from, bool stay_from,
		                      size_t base_to, bool stay_to,
			              size_t nbase){
    assert(stay_to ^ (base_from != base_to));
    return base_to * 2 * nbase + base_from + (stay_from ? nbase : 0);
}


/**  Decoding of CRF runlength model with multiple stay states
 *
 *   Runlength models have the restriction that you can't move to the same
 *   base that you are already in. Once in a move state, the model must
 *   immediately exit into either the corresponding stay state or a new
 *   move state.
 *
 *   Since this model is globally normalised, sum(exp(weights)) may not sum
 *   to one.
 *
 *   Order of params:  Each column contained 16 entries in the following order:
 *       0 --  3 : Shape parameter for discrete Weibull distribution [ACGT]
 *       4 --  7 : Scale parameter for discrete Weibull distribution [ACGT]
 *       8 -- 39 : Transition parameters
 *       ...  8 -- 23 : 4 x 4 matrix describing transition from a base-state to
 *                     base-state.  Diagonal elements are base-state to stay-state
 *       ... 24 -- 39 : 4 x 4 matrix describing transition from a stay-state to
 *                     base-state.  Diagonal elements are stay-state to stay-state
 *
 *   @param param  Flappie matrix [40 x nblk] containing predicted parameters
 *   @param path[out] Array to write out best path
 *
 *   @returns score of best path
 **/
float decode_crf_runlength(const_flappie_matrix param, int * path){
    RETURN_NULL_IF(NULL == param, NAN);
    RETURN_NULL_IF(NULL == path, NAN);
    const size_t nblk = param->nc;
    const size_t nparam = param->nr;
    const size_t nbase = nbase_from_crf_runlength_nparam(nparam);
    const size_t nstate = nbase + nbase;

    float * mem = calloc(2 * nstate, sizeof(float));
    char * traceback = calloc(nstate * nblk, sizeof(char));
    if(NULL == mem || NULL == traceback){
        free(mem);
        free(traceback);
        return NAN;
    }

    float * prev = mem;
    float * curr = mem + nstate;

    for(size_t blk=0 ; blk < nblk ; blk++){
        const size_t poffset = blk * param->stride + nbase + nbase;
        const size_t toffset = blk * nstate;
        {   //  Swap memory
            float * tmp;
            tmp = prev;
            prev = curr;
            curr = tmp;
        }

	for(size_t st=0 ; st < nstate ; st++){
            curr[st] = -HUGE_VAL;
	}

	for(size_t b1=0 ; b1 < nbase ; b1++){
            // Move into new base state.
	    for(size_t b2=0 ; b2 < nbase ; b2++){
                if(b1 == b2){
		    // Can't transition to same base
		    continue;
		}
		// From Move
		const float move_score = prev[b2] + param->data.f[poffset + rle_trans_lookup(b2, false, b1, false, nbase)];
		if(move_score > curr[b1]){
			curr[b1] = move_score;
			traceback[toffset + b1] = b2;
		}
		// From Stay
		const float stay_score = prev[b2 + nbase] + param->data.f[poffset + rle_trans_lookup(b2, true, b1, false, nbase)];
		if(stay_score > curr[b1]){
			curr[b1] = stay_score;
			traceback[toffset + b1] = b2 + nbase;
		}
	    }
	}
        for(int b=0 ; b < nbase ; b++){
            //  Stay state
            const float stay_score = prev[b + nbase] +  param->data.f[poffset + rle_trans_lookup(b, true, b, true, nbase)];
            const float move_score = prev[b] +  param->data.f[poffset + rle_trans_lookup(b, false, b, true, nbase)];
            if(stay_score > move_score){
                // Come from stay
                curr[b + nbase] = stay_score;
                traceback[toffset + b + nbase] = b + nbase;
            } else {
		// Come from move
                curr[b + nbase] = move_score;
                traceback[toffset + b + nbase] = b;
	    }
        }
    }


    size_t last_state = argmaxf(curr, nstate);
    const float logscore = curr[last_state];
    for(size_t blk=nblk ; blk > 0 ; blk--){
        const size_t blkm1 = blk - 1;
        const char state = traceback[blkm1 * nstate + last_state];
        path[blkm1] = last_state;
        last_state = state;
    }

    free(traceback);
    free(mem);

    return logscore;
}


 /**  Posterior CRF runlength model with multiple stay states
  *
  *   Runlength models have the restriction that you can't move to the same
  *   base that you are already in. Once in a move state, the model must
  *   immediately exit into either the corresponding stay state or a new
  *   move state.
  *
  *   Since this model is globally normalised, sum(exp(weights)) may not sum
  *   to one.
  *
  *   Order of params:  Each column contained 16 entries in the following order:
  *       0 --  3 : Shape parameter for discrete Weibull distribution [ACGT]
  *       4 --  7 : Scale parameter for discrete Weibull distribution [ACGT]
  *       8 -- 39 : Transition parameters
  *       ...  8 -- 23 : 4 x 4 matrix describing transition from a base-state to
  *                     base-state.  Diagonal elements are base-state to stay-state
  *       ... 24 -- 39 : 4 x 4 matrix describing transition from a stay-state to
  *                     base-state.  Diagonal elements are stay-state to stay-state
  *
  *   @param param  Flappie matrix [40 x nblk] containing predicted parameters
  *
  *   @returns score of best path
  **/
flappie_matrix transpost_crf_runlength(const_flappie_matrix param){
    RETURN_NULL_IF(NULL == param, NULL);
    const size_t nblk = param->nc;
    const size_t nparam = param->nr;
    const size_t nbase = nbase_from_crf_runlength_nparam(nparam);
    const size_t nstate = nbase + nbase;
    const size_t param_offset = nbase + nbase;

    flappie_matrix fwd = make_flappie_matrix(nstate, nblk + 1);
    flappie_matrix post = make_flappie_matrix(nparam, nblk);
    float * mem = calloc(2 * nstate, sizeof(float));
    if(NULL == fwd || NULL == post || NULL == mem){
        fwd = free_flappie_matrix(fwd);
        post = free_flappie_matrix(post);
        free(mem);
        return NULL;
    }


    for(size_t blk=0 ; blk < nblk ; blk++){
        //  Forwards calculation
        const size_t offset = blk * param->stride + param_offset;
        float * prev = fwd->data.f + blk * fwd->stride;
        float * curr = prev + fwd->stride;

        for(size_t b1=0 ; b1 < nbase ; b1++){
            curr[b1] = -HUGE_VAL;
            // Move from different base or stay
            for(size_t b2=0 ; b2 < nbase ; b2++){
                if(b1 == b2){
		    continue;
		}
		const float stay_score = prev[b2 + nbase] + param->data.f[offset + rle_trans_lookup(b2, true, b1, false, nbase)];
		const float move_score = prev[b2] + param->data.f[offset + rle_trans_lookup(b2, false, b1, false, nbase)];
		const float sum_score = logsumexpf(stay_score, move_score);
                curr[b1] = logsumexpf(curr[b1], sum_score);
            }
        }
        for(size_t b=0 ; b < nbase ; b++){
            // Stay in same base
	    const float stay_score = prev[b + nbase] + param->data.f[offset + rle_trans_lookup(b, true, b, true, nbase)];
	    const float move_score = prev[b] + param->data.f[offset + rle_trans_lookup(b, false, b, true, nbase)];
	    curr[b + nbase] = logsumexpf(stay_score, move_score);
        }
    }


    float * prev = mem;
    float * curr = prev + nstate;
    for(size_t blk=nblk ; blk > 0 ; blk--){
        const size_t foffset = (blk - 1) * fwd->stride;
        const size_t offset = (blk - 1) * param->stride + param_offset;

        // Backwards
        {  // swap
            float * tmp = curr;
            curr = prev;
            prev = tmp;
        }

        for(size_t b1=0 ; b1 < nbase ; b1++){
            curr[b1] = -HUGE_VAL;
            curr[b1 + nbase] = -HUGE_VAL;
            // Move to different base
            for(size_t b2=0 ; b2 < nbase ; b2++){
                if(b1 == b2){
		    continue;
		}
		// b1 is a move state
		const size_t move_idx = rle_trans_lookup(b1, false, b2, false, nbase);
                curr[b1] = logsumexpf(curr[b1], prev[b2] + param->data.f[offset +  move_idx]);
                post->data.f[offset + move_idx] = fwd->data.f[foffset + b1] + prev[b2] + param->data.f[offset + move_idx];
		// b1 is a stay state
		const size_t stay_idx = rle_trans_lookup(b1, true, b2, false, nbase);
                curr[b1 + nbase] = logsumexpf(curr[b1 + nbase], prev[b2] + param->data.f[offset +  stay_idx]);
                post->data.f[offset + stay_idx] = fwd->data.f[foffset + b1 + nbase] + prev[b2] + param->data.f[offset + stay_idx];
            }
        }

        for(size_t b=0 ; b < nbase ; b++){
            // Stay in same base, initial state is move
            const size_t idx = rle_trans_lookup(b, false, b, true, nbase);
            curr[b] = logsumexpf(curr[b], prev[b + nbase] + param->data.f[offset + idx]);
            post->data.f[offset + idx] = fwd->data.f[foffset + b] + param->data.f[offset + idx] + prev[b + nbase];
	}
	for(size_t b=0 ; b < nbase ; b++){
            // Stay in same base, initial state is stay
            const size_t idx = rle_trans_lookup(b, true, b, true, nbase);
            curr[b + nbase] = logsumexpf(curr[b + nbase], prev[b + nbase] + param->data.f[offset + idx]);
            post->data.f[offset + idx] = fwd->data.f[foffset + b + nbase] + param->data.f[offset + idx] + prev[b + nbase];
        }

	/*
        float score = curr[0] + fwd->data.f[foffset + 0];
        for(size_t i=1 ; i < nstate ; i++){
            score = logsumexpf(score, curr[i] + fwd->data.f[foffset + i]);
        }
        printf("score for blk %zu : %f\n", blk - 1, score);
	*/
    }

    /*
    const size_t last_offset = nblk * fwd->stride;
    float scoreF = fwd->data.f[last_offset];
    float scoreB = curr[0];
    for(size_t st=1 ; st < nstate ; st ++){
        scoreF = logsumexpf(scoreF, fwd->data.f[last_offset + st]);
        scoreB = logsumexpf(scoreB, curr[st]);
    }
    printf("ScoreF = %f  scoreB = %f\n", scoreF, scoreB);
    */

    free(mem);
    fwd = free_flappie_matrix(fwd);

    return post;
}
