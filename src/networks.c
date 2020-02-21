/*  Copyright 2018 Oxford Nanopore Technologies, Ltd */

/*  This Source Code Form is subject to the terms of the Oxford Nanopore
 *  Technologies, Ltd. Public License, v. 1.0. If a copy of the License
 *  was not  distributed with this file, You can obtain one at
 *  http://nanoporetech.com
 */

#include "layers.h"
#include "models/flipflop5_r941native.h"
#include "models/flipflop_r941native5mC.h"
#include "models/flipflop5_r103native.h"
#include "models/runlength5_r941native.h"
#include "networks.h"
#include "nnfeatures.h"
#include "flappie_stdlib.h"
#include "util.h"


enum model_type get_flappie_model_type(const char *modelstr){
    assert(NULL != modelstr);
    if(0 == strcmp(modelstr, "r941_native")) {
        return FLAPPIE_MODEL_R941_NATIVE;
    }
    if(0 == strcmp(modelstr, "r941_5mC")) {
        return FLAPPIE_MODEL_R941_5mC;
    }
    if(0 == strcmp(modelstr, "r103_native")) {
        return FLAPPIE_MODEL_R103_NATIVE;
    }
    if(0 == strcmp(modelstr, "rle_r941_native")) {
        return RUNNIE_MODEL_R941_NATIVE;
    }
    return FLAPPIE_MODEL_INVALID;
}


const char *flappie_model_string(const enum model_type model){
    switch(model){
    case FLAPPIE_MODEL_R941_NATIVE:
        return "r941_native";
    case FLAPPIE_MODEL_R941_5mC:
        return "r941_5mC";
    case FLAPPIE_MODEL_R103_NATIVE:
        return "r103_native";
    case RUNNIE_MODEL_R941_NATIVE:
        return "rle_r941_native";
    case FLAPPIE_MODEL_INVALID:
    case RUNNIE_MODEL_INVALID:
        errx(EXIT_FAILURE, "Invalid model  %s:%d", __FILE__, __LINE__);
    default:
        errx(EXIT_FAILURE, "Flappie enum failure -- report as bug. %s:%d \n", __FILE__, __LINE__);
    }
    return NULL;
}


const char *flappie_model_description(const enum model_type model){
    switch(model){
    case FLAPPIE_MODEL_R941_NATIVE:
        return "R9.4.1 model for MinION.  Trained from native DNA library";
    case FLAPPIE_MODEL_R941_5mC:
        return "R9.4.1 model for PromethION; 5mC aware.  Trained from native NA12878 library";
    case FLAPPIE_MODEL_R103_NATIVE:
        return "R10.3 model for MinION.  Trained from native DNA library";
    case RUNNIE_MODEL_R941_NATIVE:
        return "R9.4.1 run-length encoded model for MinION.  Trained from native DNA library";
    case FLAPPIE_MODEL_INVALID:
    case RUNNIE_MODEL_INVALID:
        errx(EXIT_FAILURE, "Invalid Flappie model  %s:%d", __FILE__, __LINE__);
    default:
        errx(EXIT_FAILURE, "Flappie enum failure -- report as bug. %s:%d \n", __FILE__, __LINE__);
    }
    return NULL;
}


transition_function_ptr get_transition_function(const enum model_type model){
    switch(model){
    case FLAPPIE_MODEL_R941_NATIVE:
        return flipflop5_transitions_r941native;
    case FLAPPIE_MODEL_R941_5mC:
        return flipflop_transitions_r941native5mC;
    case FLAPPIE_MODEL_R103_NATIVE:
        return flipflop5_transitions_r103native;
    case RUNNIE_MODEL_R941_NATIVE:
        return runlength5_transitions_r941native;
    case FLAPPIE_MODEL_INVALID:
    case RUNNIE_MODEL_INVALID:
        errx(EXIT_FAILURE, "Invalid Flappie model  %s:%d", __FILE__, __LINE__);
    default:
        errx(EXIT_FAILURE, "Flappie enum failure -- report as bug. %s:%d \n", __FILE__, __LINE__);
    }
    return NULL;
}


flappie_matrix calculate_transitions(const raw_table signal, float temperature, enum model_type model){
    transition_function_ptr transfun = get_transition_function(model);
    return transfun(signal, temperature);
}


typedef struct {
    //  Convolution layer
    const flappie_matrix conv_W;
    const flappie_matrix conv_b;
    int conv_stride;
    //  First modified GRU (backward)
    const flappie_matrix gruB1_iW;
    const flappie_matrix gruB1_sW;
    const flappie_matrix gruB1_sW2;
    const flappie_matrix gruB1_b;
    //  Second modified GRU (forward)
    const flappie_matrix gruF2_iW;
    const flappie_matrix gruF2_sW;
    const flappie_matrix gruF2_sW2;
    const flappie_matrix gruF2_b;
    //  Third modified GRU (backward)
    const flappie_matrix gruB3_iW;
    const flappie_matrix gruB3_sW;
    const flappie_matrix gruB3_sW2;
    const flappie_matrix gruB3_b;
    //  Fourth modified GRU (forward)
    const flappie_matrix gruF4_iW;
    const flappie_matrix gruF4_sW;
    const flappie_matrix gruF4_sW2;
    const flappie_matrix gruF4_b;
    //  Fifth modified GRU (backward)
    const flappie_matrix gruB5_iW;
    const flappie_matrix gruB5_sW;
    const flappie_matrix gruB5_sW2;
    const flappie_matrix gruB5_b;
    //  Output
    const flappie_matrix FF_W;
    const flappie_matrix FF_b;
} sloika_model;


typedef struct {
    //  Convolution layer
    const flappie_matrix conv_W;
    const flappie_matrix conv_b;
    int conv_stride;
    //  First modified GRU (backward)
    const flappie_matrix gruB1_iW;
    const flappie_matrix gruB1_sW;
    const flappie_matrix gruB1_b;
    //  Second modified GRU (forward)
    const flappie_matrix gruF2_iW;
    const flappie_matrix gruF2_sW;
    const flappie_matrix gruF2_b;
    //  Third modified GRU (backward)
    const flappie_matrix gruB3_iW;
    const flappie_matrix gruB3_sW;
    const flappie_matrix gruB3_b;
    //  Fourth modified GRU (forward)
    const flappie_matrix gruF4_iW;
    const flappie_matrix gruF4_sW;
    const flappie_matrix gruF4_b;
    //  Fifth modified GRU (backward)
    const flappie_matrix gruB5_iW;
    const flappie_matrix gruB5_sW;
    const flappie_matrix gruB5_b;
    //  Output
    const flappie_matrix FF_W;
    const flappie_matrix FF_b;
} guppy_model;


typedef struct {
    //  Convolution layer
    const flappie_matrix conv1_W;
    const flappie_matrix conv1_b;
    int conv1_stride;
    const flappie_matrix conv2_W;
    const flappie_matrix conv2_b;
    int conv2_stride;
    const flappie_matrix conv3_W;
    const flappie_matrix conv3_b;
    int conv3_stride;
    //  First modified LSTM (backward)
    const flappie_matrix lstmB1_iW;
    const flappie_matrix lstmB1_sW;
    const flappie_matrix lstmB1_b;
    //  Second modified LSTM (forward)
    const flappie_matrix lstmF2_iW;
    const flappie_matrix lstmF2_sW;
    const flappie_matrix lstmF2_b;
    //  Third modified LSTM (backward)
    const flappie_matrix lstmB3_iW;
    const flappie_matrix lstmB3_sW;
    const flappie_matrix lstmB3_b;
    //  Fourth modified LSTM (forward)
    const flappie_matrix lstmF4_iW;
    const flappie_matrix lstmF4_sW;
    const flappie_matrix lstmF4_b;
    //  Fifth modified LSTM (backward)
    const flappie_matrix lstmB5_iW;
    const flappie_matrix lstmB5_sW;
    const flappie_matrix lstmB5_b;
    //  Output
    const flappie_matrix FF_W;
    const flappie_matrix FF_b;
} guppy_stride5_model;


guppy_stride5_model flipflop5_r941native_guppy = {
    //  Convolution layer
    .conv1_W = &_conv1_rnnrf_flipflop5_r941native_W,
    .conv1_b = &_conv1_rnnrf_flipflop5_r941native_b,
    .conv1_stride = conv1_rnnrf_flipflop5_r941native_stride,
    .conv2_W = &_conv2_rnnrf_flipflop5_r941native_W,
    .conv2_b = &_conv2_rnnrf_flipflop5_r941native_b,
    .conv2_stride = conv2_rnnrf_flipflop5_r941native_stride,
    .conv3_W = &_conv3_rnnrf_flipflop5_r941native_W,
    .conv3_b = &_conv3_rnnrf_flipflop5_r941native_b,
    .conv3_stride = conv3_rnnrf_flipflop5_r941native_stride,
    //.conv_stride = 2,
    //  First modified LSTM (backward)
    .lstmB1_iW = &_lstmB1_rnnrf_flipflop5_r941native_iW,
    .lstmB1_sW = &_lstmB1_rnnrf_flipflop5_r941native_sW,
    .lstmB1_b = &_lstmB1_rnnrf_flipflop5_r941native_b,
    //  Second modified LSTM (forward)
    .lstmF2_iW = &_lstmF2_rnnrf_flipflop5_r941native_iW,
    .lstmF2_sW = &_lstmF2_rnnrf_flipflop5_r941native_sW,
    .lstmF2_b = &_lstmF2_rnnrf_flipflop5_r941native_b,
    //  Third modified LSTM (backward)
    .lstmB3_iW = &_lstmB3_rnnrf_flipflop5_r941native_iW,
    .lstmB3_sW = &_lstmB3_rnnrf_flipflop5_r941native_sW,
    .lstmB3_b = &_lstmB3_rnnrf_flipflop5_r941native_b,
    //  Fourth modified LSTM (forward)
    .lstmF4_iW = &_lstmF4_rnnrf_flipflop5_r941native_iW,
    .lstmF4_sW = &_lstmF4_rnnrf_flipflop5_r941native_sW,
    .lstmF4_b = &_lstmF4_rnnrf_flipflop5_r941native_b,
    //  Fifth modified LSTM (backward)
    .lstmB5_iW = &_lstmB5_rnnrf_flipflop5_r941native_iW,
    .lstmB5_sW = &_lstmB5_rnnrf_flipflop5_r941native_sW,
    .lstmB5_b = &_lstmB5_rnnrf_flipflop5_r941native_b,
    //  Output
    .FF_W = &_FF_rnnrf_flipflop5_r941native_W,
    .FF_b = &_FF_rnnrf_flipflop5_r941native_b
};


guppy_model flipflop_r941native5mC_guppy = {
    //  Convolution layer
    .conv_W = &_conv_rnnrf_flipflop_r941native5mC_W,
    .conv_b = &_conv_rnnrf_flipflop_r941native5mC_b,
    .conv_stride = conv_rnnrf_flipflop_r941native5mC_stride,
    //.conv_stride = 2,
    //  First modified GRU (backward)
    .gruB1_iW = &_gruB1_rnnrf_flipflop_r941native5mC_iW,
    .gruB1_sW = &_gruB1_rnnrf_flipflop_r941native5mC_sW,
    .gruB1_b = &_gruB1_rnnrf_flipflop_r941native5mC_b,
    //  Second modified GRU (forward)
    .gruF2_iW = &_gruF2_rnnrf_flipflop_r941native5mC_iW,
    .gruF2_sW = &_gruF2_rnnrf_flipflop_r941native5mC_sW,
    .gruF2_b = &_gruF2_rnnrf_flipflop_r941native5mC_b,
    //  Third modified GRU (backward)
    .gruB3_iW = &_gruB3_rnnrf_flipflop_r941native5mC_iW,
    .gruB3_sW = &_gruB3_rnnrf_flipflop_r941native5mC_sW,
    .gruB3_b = &_gruB3_rnnrf_flipflop_r941native5mC_b,
    //  Fourth modified GRU (forward)
    .gruF4_iW = &_gruF4_rnnrf_flipflop_r941native5mC_iW,
    .gruF4_sW = &_gruF4_rnnrf_flipflop_r941native5mC_sW,
    .gruF4_b = &_gruF4_rnnrf_flipflop_r941native5mC_b,
    //  Fifth modified GRU (backward)
    .gruB5_iW = &_gruB5_rnnrf_flipflop_r941native5mC_iW,
    .gruB5_sW = &_gruB5_rnnrf_flipflop_r941native5mC_sW,
    .gruB5_b = &_gruB5_rnnrf_flipflop_r941native5mC_b,
    //  Output
    .FF_W = &_FF_rnnrf_flipflop_r941native5mC_W,
    .FF_b = &_FF_rnnrf_flipflop_r941native5mC_b
};


guppy_stride5_model flipflop5_r103native_guppy = {
    //  Convolution layer
    .conv1_W = &_conv1_rnnrf_flipflop5_r103native_W,
    .conv1_b = &_conv1_rnnrf_flipflop5_r103native_b,
    .conv1_stride = conv1_rnnrf_flipflop5_r103native_stride,
    .conv2_W = &_conv2_rnnrf_flipflop5_r103native_W,
    .conv2_b = &_conv2_rnnrf_flipflop5_r103native_b,
    .conv2_stride = conv2_rnnrf_flipflop5_r103native_stride,
    .conv3_W = &_conv3_rnnrf_flipflop5_r103native_W,
    .conv3_b = &_conv3_rnnrf_flipflop5_r103native_b,
    .conv3_stride = conv3_rnnrf_flipflop5_r103native_stride,
    //.conv_stride = 2,
    //  First modified LSTM (backward)
    .lstmB1_iW = &_lstmB1_rnnrf_flipflop5_r103native_iW,
    .lstmB1_sW = &_lstmB1_rnnrf_flipflop5_r103native_sW,
    .lstmB1_b = &_lstmB1_rnnrf_flipflop5_r103native_b,
    //  Second modified LSTM (forward)
    .lstmF2_iW = &_lstmF2_rnnrf_flipflop5_r103native_iW,
    .lstmF2_sW = &_lstmF2_rnnrf_flipflop5_r103native_sW,
    .lstmF2_b = &_lstmF2_rnnrf_flipflop5_r103native_b,
    //  Third modified LSTM (backward)
    .lstmB3_iW = &_lstmB3_rnnrf_flipflop5_r103native_iW,
    .lstmB3_sW = &_lstmB3_rnnrf_flipflop5_r103native_sW,
    .lstmB3_b = &_lstmB3_rnnrf_flipflop5_r103native_b,
    //  Fourth modified LSTM (forward)
    .lstmF4_iW = &_lstmF4_rnnrf_flipflop5_r103native_iW,
    .lstmF4_sW = &_lstmF4_rnnrf_flipflop5_r103native_sW,
    .lstmF4_b = &_lstmF4_rnnrf_flipflop5_r103native_b,
    //  Fifth modified LSTM (backward)
    .lstmB5_iW = &_lstmB5_rnnrf_flipflop5_r103native_iW,
    .lstmB5_sW = &_lstmB5_rnnrf_flipflop5_r103native_sW,
    .lstmB5_b = &_lstmB5_rnnrf_flipflop5_r103native_b,
    //  Output
    .FF_W = &_FF_rnnrf_flipflop5_r103native_W,
    .FF_b = &_FF_rnnrf_flipflop5_r103native_b
};


guppy_stride5_model runlength5_r941native_guppy = {
    //  Convolution layer
    .conv1_W = &_conv1_rnnrf_rle5_r941native_W,
    .conv1_b = &_conv1_rnnrf_rle5_r941native_b,
    .conv1_stride = conv1_rnnrf_rle5_r941native_stride,
    .conv2_W = &_conv2_rnnrf_rle5_r941native_W,
    .conv2_b = &_conv2_rnnrf_rle5_r941native_b,
    .conv2_stride = conv2_rnnrf_rle5_r941native_stride,
    .conv3_W = &_conv3_rnnrf_rle5_r941native_W,
    .conv3_b = &_conv3_rnnrf_rle5_r941native_b,
    .conv3_stride = conv3_rnnrf_rle5_r941native_stride,
    //.conv_stride = 2,
    //  First modified LSTM (backward)
    .lstmB1_iW = &_lstmB1_rnnrf_rle5_r941native_iW,
    .lstmB1_sW = &_lstmB1_rnnrf_rle5_r941native_sW,
    .lstmB1_b = &_lstmB1_rnnrf_rle5_r941native_b,
    //  Second modified LSTM (forward)
    .lstmF2_iW = &_lstmF2_rnnrf_rle5_r941native_iW,
    .lstmF2_sW = &_lstmF2_rnnrf_rle5_r941native_sW,
    .lstmF2_b = &_lstmF2_rnnrf_rle5_r941native_b,
    //  Third modified LSTM (backward)
    .lstmB3_iW = &_lstmB3_rnnrf_rle5_r941native_iW,
    .lstmB3_sW = &_lstmB3_rnnrf_rle5_r941native_sW,
    .lstmB3_b = &_lstmB3_rnnrf_rle5_r941native_b,
    //  Fourth modified LSTM (forward)
    .lstmF4_iW = &_lstmF4_rnnrf_rle5_r941native_iW,
    .lstmF4_sW = &_lstmF4_rnnrf_rle5_r941native_sW,
    .lstmF4_b = &_lstmF4_rnnrf_rle5_r941native_b,
    //  Fifth modified LSTM (backward)
    .lstmB5_iW = &_lstmB5_rnnrf_rle5_r941native_iW,
    .lstmB5_sW = &_lstmB5_rnnrf_rle5_r941native_sW,
    .lstmB5_b = &_lstmB5_rnnrf_rle5_r941native_b,
    //  Output
    .FF_W = &_FF_rnnrf_rle5_r941native_W,
    .FF_b = &_FF_rnnrf_rle5_r941native_b
};



flappie_matrix flipflop_gru_transitions(const raw_table signal, float temperature, const sloika_model * net){
    RETURN_NULL_IF(0 == signal.n, NULL);
    RETURN_NULL_IF(NULL == signal.raw, NULL);

    flappie_matrix raw_mat = features_from_raw(signal);
    flappie_matrix conv =
        convolution(raw_mat, net->conv_W, net->conv_b, net->conv_stride, NULL);
    elu_activation_inplace(conv);
    raw_mat = free_flappie_matrix(raw_mat);
    //  First GRU layer
    flappie_matrix gruB1in = feedforward_linear(conv, net->gruB1_iW, net->gruB1_b, NULL);
    flappie_matrix gruB1 = gru_backward(gruB1in, net->gruB1_sW, net->gruB1_sW2, NULL);
    residual_inplace(conv, gruB1);
    conv = free_flappie_matrix(conv);
    gruB1in = free_flappie_matrix(gruB1in);
    //  Second GRU layer
    flappie_matrix gruF2in = feedforward_linear(gruB1, net->gruF2_iW, net->gruF2_b, NULL);
    flappie_matrix gruF2 = gru_forward(gruF2in, net->gruF2_sW, net->gruF2_sW2, NULL);
    residual_inplace(gruB1, gruF2);
    gruB1 = free_flappie_matrix(gruB1);
    gruF2in = free_flappie_matrix(gruF2in);
    //  Third GRU layer
    flappie_matrix gruB3in = feedforward_linear(gruF2, net->gruB3_iW, net->gruB3_b, NULL);
    flappie_matrix gruB3 = gru_backward(gruB3in, net->gruB3_sW, net->gruB3_sW2, NULL);
    residual_inplace(gruF2, gruB3);
    gruF2 = free_flappie_matrix(gruF2);
    gruB3in = free_flappie_matrix(gruB3in);
    //  Fourth GRU layer
    flappie_matrix gruF4in = feedforward_linear(gruB3, net->gruF4_iW, net->gruF4_b, NULL);
    flappie_matrix gruF4 = gru_forward(gruF4in, net->gruF4_sW, net->gruF4_sW2, NULL);
    residual_inplace(gruB3, gruF4);
    gruB3 = free_flappie_matrix(gruB3);
    gruF4in = free_flappie_matrix(gruF4in);
    //  Fifth GRU layer
    flappie_matrix gruB5in = feedforward_linear(gruF4, net->gruB5_iW, net->gruB5_b, NULL);
    flappie_matrix gruB5 = gru_backward(gruB5in, net->gruB5_sW, net->gruB5_sW2, NULL);
    residual_inplace(gruF4, gruB5);
    gruF4 = free_flappie_matrix(gruF4);
    gruB5in = free_flappie_matrix(gruB5in);

    flappie_matrix trans = globalnorm_flipflop(gruB5, net->FF_W, net->FF_b, temperature, NULL);
    gruB5 = free_flappie_matrix(gruB5);

    return trans;
}


flappie_matrix flipflop_guppy_transitions(const raw_table signal, float temperature, const guppy_model * net){
    RETURN_NULL_IF(0 == signal.n, NULL);
    RETURN_NULL_IF(NULL == signal.raw, NULL);

    flappie_matrix raw_mat = features_from_raw(signal);
    flappie_matrix conv =
        convolution(raw_mat, net->conv_W, net->conv_b, net->conv_stride, NULL);
    tanh_activation_inplace(conv);
    raw_mat = free_flappie_matrix(raw_mat);
    //  First GRU layer
    flappie_matrix gruB1in = feedforward_linear(conv, net->gruB1_iW, net->gruB1_b, NULL);
    conv = free_flappie_matrix(conv);
    flappie_matrix gruB1 = grumod_backward(gruB1in, net->gruB1_sW, NULL);
    gruB1in = free_flappie_matrix(gruB1in);
    //  Second GRU layer
    flappie_matrix gruF2in = feedforward_linear(gruB1, net->gruF2_iW, net->gruF2_b, NULL);
    gruB1 = free_flappie_matrix(gruB1);
    flappie_matrix gruF2 = grumod_forward(gruF2in, net->gruF2_sW, NULL);
    gruF2in = free_flappie_matrix(gruF2in);
    //  Third GRU layer
    flappie_matrix gruB3in = feedforward_linear(gruF2, net->gruB3_iW, net->gruB3_b, NULL);
    gruF2 = free_flappie_matrix(gruF2);
    flappie_matrix gruB3 = grumod_backward(gruB3in, net->gruB3_sW, NULL);
    gruB3in = free_flappie_matrix(gruB3in);
    //  Fourth GRU layer
    flappie_matrix gruF4in = feedforward_linear(gruB3, net->gruF4_iW, net->gruF4_b, NULL);
    gruB3 = free_flappie_matrix(gruB3);
    flappie_matrix gruF4 = grumod_forward(gruF4in, net->gruF4_sW, NULL);
    gruF4in = free_flappie_matrix(gruF4in);
    //  Fifth GRU layer
    flappie_matrix gruB5in = feedforward_linear(gruF4, net->gruB5_iW, net->gruB5_b, NULL);
    gruF4 = free_flappie_matrix(gruF4);
    flappie_matrix gruB5 = grumod_backward(gruB5in, net->gruB5_sW, NULL);
    gruB5in = free_flappie_matrix(gruB5in);

    flappie_matrix trans = globalnorm_flipflop(gruB5, net->FF_W, net->FF_b, temperature, NULL);
    gruB5 = free_flappie_matrix(gruB5);

    return trans;
}


flappie_matrix flipflop_relu_transitions(const raw_table signal, float temperature, const sloika_model * net){
    RETURN_NULL_IF(0 == signal.n, NULL);
    RETURN_NULL_IF(NULL == signal.raw, NULL);

    flappie_matrix raw_mat = features_from_raw(signal);
    flappie_matrix conv =
        convolution(raw_mat, net->conv_W, net->conv_b, net->conv_stride, NULL);
    elu_activation_inplace(conv);
    raw_mat = free_flappie_matrix(raw_mat);
    //  First GRU layer
    flappie_matrix gruB1in = feedforward_linear(conv, net->gruB1_iW, net->gruB1_b, NULL);
    flappie_matrix gruB1 = gru_relu_backward(gruB1in, net->gruB1_sW, net->gruB1_sW2, NULL);
    residual_inplace(conv, gruB1);
    conv = free_flappie_matrix(conv);
    gruB1in = free_flappie_matrix(gruB1in);
    //  Second GRU layer
    flappie_matrix gruF2in = feedforward_linear(gruB1, net->gruF2_iW, net->gruF2_b, NULL);
    flappie_matrix gruF2 = gru_relu_forward(gruF2in, net->gruF2_sW, net->gruF2_sW2, NULL);
    residual_inplace(gruB1, gruF2);
    gruB1 = free_flappie_matrix(gruB1);
    gruF2in = free_flappie_matrix(gruF2in);
    //  Third GRU layer
    flappie_matrix gruB3in = feedforward_linear(gruF2, net->gruB3_iW, net->gruB3_b, NULL);
    flappie_matrix gruB3 = gru_relu_backward(gruB3in, net->gruB3_sW, net->gruB3_sW2, NULL);
    residual_inplace(gruF2, gruB3);
    gruF2 = free_flappie_matrix(gruF2);
    gruB3in = free_flappie_matrix(gruB3in);
    //  Fourth GRU layer
    flappie_matrix gruF4in = feedforward_linear(gruB3, net->gruF4_iW, net->gruF4_b, NULL);
    flappie_matrix gruF4 = gru_relu_forward(gruF4in, net->gruF4_sW, net->gruF4_sW2, NULL);
    residual_inplace(gruB3, gruF4);
    gruB3 = free_flappie_matrix(gruB3);
    gruF4in = free_flappie_matrix(gruF4in);
    //  Fifth GRU layer
    flappie_matrix gruB5in = feedforward_linear(gruF4, net->gruB5_iW, net->gruB5_b, NULL);
    flappie_matrix gruB5 = gru_relu_backward(gruB5in, net->gruB5_sW, net->gruB5_sW2, NULL);
    residual_inplace(gruF4, gruB5);
    gruF4 = free_flappie_matrix(gruF4);
    gruB5in = free_flappie_matrix(gruB5in);

    flappie_matrix trans = globalnorm_flipflop(gruB5, net->FF_W, net->FF_b, temperature, NULL);
    gruB5 = free_flappie_matrix(gruB5);

    return trans;
}


flappie_matrix flipflop5_guppy_transitions(const raw_table signal, float temperature, const guppy_stride5_model * net){
    RETURN_NULL_IF(0 == signal.n, NULL);
    RETURN_NULL_IF(NULL == signal.raw, NULL);

    flappie_matrix raw_mat = features_from_raw(signal);
    flappie_matrix conv1 =
        convolution(raw_mat, net->conv1_W, net->conv1_b, net->conv1_stride, NULL);
    swish_activation_inplace(conv1);
    raw_mat = free_flappie_matrix(raw_mat);
    flappie_matrix conv2 =
        convolution(conv1, net->conv2_W, net->conv2_b, net->conv2_stride, NULL);
    swish_activation_inplace(conv2);
    conv1 = free_flappie_matrix(conv1);
    flappie_matrix conv3 =
        convolution(conv2, net->conv3_W, net->conv3_b, net->conv3_stride, NULL);
    swish_activation_inplace(conv3);
    conv2 = free_flappie_matrix(conv2);
    //  First LSTM layer
    flappie_matrix lstmB1in = feedforward_linear(conv3, net->lstmB1_iW, net->lstmB1_b, NULL);
    conv3 = free_flappie_matrix(conv3);
    flappie_matrix lstmB1 = lstm_backward(lstmB1in, net->lstmB1_sW, NULL);
    lstmB1in = free_flappie_matrix(lstmB1in);
    //  Second LSTM layer
    flappie_matrix lstmF2in = feedforward_linear(lstmB1, net->lstmF2_iW, net->lstmF2_b, NULL);
    lstmB1 = free_flappie_matrix(lstmB1);
    flappie_matrix lstmF2 = lstm_forward(lstmF2in, net->lstmF2_sW, NULL);
    lstmF2in = free_flappie_matrix(lstmF2in);
    //  Third LSTM layer
    flappie_matrix lstmB3in = feedforward_linear(lstmF2, net->lstmB3_iW, net->lstmB3_b, NULL);
    lstmF2 = free_flappie_matrix(lstmF2);
    flappie_matrix lstmB3 = lstm_backward(lstmB3in, net->lstmB3_sW, NULL);
    lstmB3in = free_flappie_matrix(lstmB3in);
    //  Fourth LSTM layer
    flappie_matrix lstmF4in = feedforward_linear(lstmB3, net->lstmF4_iW, net->lstmF4_b, NULL);
    lstmB3 = free_flappie_matrix(lstmB3);
    flappie_matrix lstmF4 = lstm_forward(lstmF4in, net->lstmF4_sW, NULL);
    lstmF4in = free_flappie_matrix(lstmF4in);
    //  Fifth LSTM layer
    flappie_matrix lstmB5in = feedforward_linear(lstmF4, net->lstmB5_iW, net->lstmB5_b, NULL);
    lstmF4 = free_flappie_matrix(lstmF4);
    flappie_matrix lstmB5 = lstm_backward(lstmB5in, net->lstmB5_sW, NULL);
    lstmB5in = free_flappie_matrix(lstmB5in);

    flappie_matrix trans = globalnorm_flipflop(lstmB5, net->FF_W, net->FF_b, temperature, NULL);
    lstmB5 = free_flappie_matrix(lstmB5);

    return trans;
}


flappie_matrix runlength_guppy_transitions(const raw_table signal, float temperature, const guppy_model * net){
    RETURN_NULL_IF(0 == signal.n, NULL);
    RETURN_NULL_IF(NULL == signal.raw, NULL);

    flappie_matrix raw_mat = features_from_raw(signal);
    flappie_matrix conv =
        convolution(raw_mat, net->conv_W, net->conv_b, net->conv_stride, NULL);
    tanh_activation_inplace(conv);
    raw_mat = free_flappie_matrix(raw_mat);
    //  First GRU layer
    flappie_matrix gruB1in = feedforward_linear(conv, net->gruB1_iW, net->gruB1_b, NULL);
    conv = free_flappie_matrix(conv);
    flappie_matrix gruB1 = grumod_backward(gruB1in, net->gruB1_sW, NULL);
    gruB1in = free_flappie_matrix(gruB1in);
    //  Second GRU layer
    flappie_matrix gruF2in = feedforward_linear(gruB1, net->gruF2_iW, net->gruF2_b, NULL);
    gruB1 = free_flappie_matrix(gruB1);
    flappie_matrix gruF2 = grumod_forward(gruF2in, net->gruF2_sW, NULL);
    gruF2in = free_flappie_matrix(gruF2in);
    //  Third GRU layer
    flappie_matrix gruB3in = feedforward_linear(gruF2, net->gruB3_iW, net->gruB3_b, NULL);
    gruF2 = free_flappie_matrix(gruF2);
    flappie_matrix gruB3 = grumod_backward(gruB3in, net->gruB3_sW, NULL);
    gruB3in = free_flappie_matrix(gruB3in);
    //  Fourth GRU layer
    flappie_matrix gruF4in = feedforward_linear(gruB3, net->gruF4_iW, net->gruF4_b, NULL);
    gruB3 = free_flappie_matrix(gruB3);
    flappie_matrix gruF4 = grumod_forward(gruF4in, net->gruF4_sW, NULL);
    gruF4in = free_flappie_matrix(gruF4in);
    //  Fifth GRU layer
    flappie_matrix gruB5in = feedforward_linear(gruF4, net->gruB5_iW, net->gruB5_b, NULL);
    gruF4 = free_flappie_matrix(gruF4);
    flappie_matrix gruB5 = grumod_backward(gruB5in, net->gruB5_sW, NULL);
    gruB5in = free_flappie_matrix(gruB5in);

    flappie_matrix trans = globalnorm_runlength(gruB5, net->FF_W, net->FF_b, temperature, NULL);

    gruB5 = free_flappie_matrix(gruB5);

    return trans;
}


flappie_matrix runlengthV2_guppy_transitions(const raw_table signal, float temperature, const guppy_model * net){
    RETURN_NULL_IF(0 == signal.n, NULL);
    RETURN_NULL_IF(NULL == signal.raw, NULL);

    flappie_matrix raw_mat = features_from_raw(signal);
    flappie_matrix conv =
        convolution(raw_mat, net->conv_W, net->conv_b, net->conv_stride, NULL);
    tanh_activation_inplace(conv);
    raw_mat = free_flappie_matrix(raw_mat);
    //  First GRU layer
    flappie_matrix gruB1in = feedforward_linear(conv, net->gruB1_iW, net->gruB1_b, NULL);
    conv = free_flappie_matrix(conv);
    flappie_matrix gruB1 = lstm_backward(gruB1in, net->gruB1_sW, NULL);
    gruB1in = free_flappie_matrix(gruB1in);
    //  Second GRU layer
    flappie_matrix gruF2in = feedforward_linear(gruB1, net->gruF2_iW, net->gruF2_b, NULL);
    gruB1 = free_flappie_matrix(gruB1);
    flappie_matrix gruF2 = lstm_forward(gruF2in, net->gruF2_sW, NULL);
    gruF2in = free_flappie_matrix(gruF2in);
    //  Third GRU layer
    flappie_matrix gruB3in = feedforward_linear(gruF2, net->gruB3_iW, net->gruB3_b, NULL);
    gruF2 = free_flappie_matrix(gruF2);
    flappie_matrix gruB3 = lstm_backward(gruB3in, net->gruB3_sW, NULL);
    gruB3in = free_flappie_matrix(gruB3in);
    //  Fourth GRU layer
    flappie_matrix gruF4in = feedforward_linear(gruB3, net->gruF4_iW, net->gruF4_b, NULL);
    gruB3 = free_flappie_matrix(gruB3);
    flappie_matrix gruF4 = lstm_forward(gruF4in, net->gruF4_sW, NULL);
    gruF4in = free_flappie_matrix(gruF4in);
    //  Fifth GRU layer
    flappie_matrix gruB5in = feedforward_linear(gruF4, net->gruB5_iW, net->gruB5_b, NULL);
    gruF4 = free_flappie_matrix(gruF4);
    flappie_matrix gruB5 = lstm_backward(gruB5in, net->gruB5_sW, NULL);
    gruB5in = free_flappie_matrix(gruB5in);

    flappie_matrix trans = globalnorm_runlengthV2(gruB5, net->FF_W, net->FF_b, temperature, NULL);

    gruB5 = free_flappie_matrix(gruB5);

    return trans;
}


flappie_matrix runlength5_guppy_transitions(const raw_table signal, float temperature, const guppy_stride5_model * net){
    RETURN_NULL_IF(0 == signal.n, NULL);
    RETURN_NULL_IF(NULL == signal.raw, NULL);

    flappie_matrix raw_mat = features_from_raw(signal);
    flappie_matrix conv1 =
        convolution(raw_mat, net->conv1_W, net->conv1_b, net->conv1_stride, NULL);
    swish_activation_inplace(conv1);
    raw_mat = free_flappie_matrix(raw_mat);
    flappie_matrix conv2 =
        convolution(conv1, net->conv2_W, net->conv2_b, net->conv2_stride, NULL);
    swish_activation_inplace(conv2);
    conv1 = free_flappie_matrix(conv1);
    flappie_matrix conv3 =
        convolution(conv2, net->conv3_W, net->conv3_b, net->conv3_stride, NULL);
    swish_activation_inplace(conv3);
    conv2 = free_flappie_matrix(conv2);
    //  First LSTM layer
    flappie_matrix lstmB1in = feedforward_linear(conv3, net->lstmB1_iW, net->lstmB1_b, NULL);
    conv3 = free_flappie_matrix(conv3);
    flappie_matrix lstmB1 = lstm_backward(lstmB1in, net->lstmB1_sW, NULL);
    lstmB1in = free_flappie_matrix(lstmB1in);
    //  Second LSTM layer
    flappie_matrix lstmF2in = feedforward_linear(lstmB1, net->lstmF2_iW, net->lstmF2_b, NULL);
    lstmB1 = free_flappie_matrix(lstmB1);
    flappie_matrix lstmF2 = lstm_forward(lstmF2in, net->lstmF2_sW, NULL);
    lstmF2in = free_flappie_matrix(lstmF2in);
    //  Third LSTM layer
    flappie_matrix lstmB3in = feedforward_linear(lstmF2, net->lstmB3_iW, net->lstmB3_b, NULL);
    lstmF2 = free_flappie_matrix(lstmF2);
    flappie_matrix lstmB3 = lstm_backward(lstmB3in, net->lstmB3_sW, NULL);
    lstmB3in = free_flappie_matrix(lstmB3in);
    //  Fourth LSTM layer
    flappie_matrix lstmF4in = feedforward_linear(lstmB3, net->lstmF4_iW, net->lstmF4_b, NULL);
    lstmB3 = free_flappie_matrix(lstmB3);
    flappie_matrix lstmF4 = lstm_forward(lstmF4in, net->lstmF4_sW, NULL);
    lstmF4in = free_flappie_matrix(lstmF4in);
    //  Fifth LSTM layer
    flappie_matrix lstmB5in = feedforward_linear(lstmF4, net->lstmB5_iW, net->lstmB5_b, NULL);
    lstmF4 = free_flappie_matrix(lstmF4);
    flappie_matrix lstmB5 = lstm_backward(lstmB5in, net->lstmB5_sW, NULL);
    lstmB5in = free_flappie_matrix(lstmB5in);

    flappie_matrix trans = globalnorm_runlengthV2(lstmB5, net->FF_W, net->FF_b, temperature, NULL);
    lstmB5 = free_flappie_matrix(lstmB5);

    return trans;
}


flappie_matrix flipflop5_transitions_r941native(const raw_table signal, float temperature){
    return flipflop5_guppy_transitions(signal, temperature, &flipflop5_r941native_guppy);
}

flappie_matrix flipflop_transitions_r941native5mC(const raw_table signal, float temperature){
    return flipflop_guppy_transitions(signal, temperature, &flipflop_r941native5mC_guppy);
}

flappie_matrix flipflop5_transitions_r103native(const raw_table signal, float temperature){
    return flipflop5_guppy_transitions(signal, temperature, &flipflop5_r103native_guppy);
}

flappie_matrix runlength5_transitions_r941native(const raw_table signal, float temperature){
    return runlength5_guppy_transitions(signal, temperature, &runlength5_r941native_guppy);
}
