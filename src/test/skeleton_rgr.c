/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. License, v. 1.0. If a copy of the License was not
 *  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#include "fast5_interface.h"
#include "layers.h"
#include "nanonet_rgr.h"
#include "nnfeatures.h"
#include "flappie_common.h"
#include "flappie_util.h"

int main(int argc, char * argv[]){
    int trim_start = 200;
    int trim_end = 10;
    int varseg_chunk = 100;
    float varseg_thresh = 0;

    char * filename = argv[1];


    raw_table signal = read_raw(filename, true);
    signal =  trim_and_segment_raw(signal, trim_start, trim_end, varseg_chunk, varseg_thresh);


    flappie_matrix raw_mat = features_from_raw(signal);
    write_flappie_matrix("input.crp", raw_mat);
    flappie_matrix conv =
        convolution(raw_mat, conv_rgr_W, conv_rgr_b, conv_rgr_stride, NULL);
    elu_activation_inplace(conv);
    write_flappie_matrix("convolution.crp", conv);


    raw_mat = free_flappie_matrix(raw_mat);

    //  First GRU layer
    flappie_matrix gruB1 =
        gru_backward(conv, gruB1_rgr_iW, gruB1_rgr_sW, gruB1_rgr_sW2,
                     gruB1_rgr_b, NULL);
    write_flappie_matrix("gruB1.crp", gruB1);
    conv = free_flappie_matrix(conv);
    //  Second GRU layer
    flappie_matrix gruF2 =
        gru_forward(gruB1, gruF2_rgr_iW, gruF2_rgr_sW, gruF2_rgr_sW2,
                    gruF2_rgr_b, NULL);
    write_flappie_matrix("gruF2.crp", gruF2);
    gruB1 = free_flappie_matrix(gruB1);
    //  Thrid GRU layer
    flappie_matrix gruB3 =
        gru_backward(gruF2, gruB3_rgr_iW, gruB3_rgr_sW, gruB3_rgr_sW2,
                     gruB3_rgr_b, NULL);
    write_flappie_matrix("gruB3.crp", gruB3);
    gruF2 = free_flappie_matrix(gruF2);

    flappie_matrix post = softmax(gruB3, FF_rgr_W, FF_rgr_b, NULL);
    write_flappie_matrix("softmax.crp", post);
    gruB3 = free_flappie_matrix(gruB3);

    post = free_flappie_matrix(post);
}
