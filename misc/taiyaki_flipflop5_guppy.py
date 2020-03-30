#!/usr/bin/env python3

#  Copyright 2018 Oxford Nanopore Technologies, Ltd

#  This Source Code Form is subject to the terms of the Oxford Nanopore
#  Technologies, Ltd. Public License, v. 1.0. If a copy of the License
#  was not  distributed with this file, You can obtain one at
#  http://nanoporetech.com

import argparse
import math
import re
import sys

from taiyaki import helpers
from taiyaki.cmdargs import AutoBool, FileExists
from taiyaki.layers import _cudnn_to_guppy_gru, DeltaSample

parser = argparse.ArgumentParser()
parser.add_argument('--id', default='' , help='Identifier for model names')
parser.add_argument('--scale', default=False, action=AutoBool,
                    help='Correct scaling when network trained without MAD factor')
parser.add_argument('model', action=FileExists, help='Pickle to read model from')


trim_trailing_zeros = re.compile('0+p')

def small_hex(f):
    hf = float(f).hex()
    return trim_trailing_zeros.sub('p', hf)


def process_column(v, pad):
    """ process and pad """
    return [small_hex(f) for f in v] + [small_hex(0.0)] * pad


def cformatM(fh, name, X, nr=None, nc=None):
    nrq = int(math.ceil(X.shape[1] / 4.0))
    pad = nrq * 4 - X.shape[1]
    lines = map(lambda v: ', '.join(process_column(v, pad)), X)

    if nr is None:
        nr = X.shape[1]
    else:
        nrq = int(math.ceil(nr / 4.0))
    if nc is None:
        nc = X.shape[0]

    fh.write('float {}[] = {}\n'.format('__' + name, '{'))
    fh.write('\t' + ',\n\t'.join(lines))
    fh.write('};\n')
    fh.write('_Mat {} = {}\n\t.nr = {},\n\t.nrq = {},\n\t.nc = {},\n\t.stride = {},\n\t.data.f = {}\n{};\n'.format('_' + name, '{', nr, nrq, nc, nrq * 4, '__' + name, '}'))
    fh.write('const flappie_matrix {} = &{};\n\n'.format(name, '_' + name))


def cformatV(fh, name, X):
    nrq = int(math.ceil(X.shape[0] / 4.0))
    pad = nrq * 4 - X.shape[0]
    lines = ', '.join(list(map(lambda f: small_hex(f), X)) + [small_hex(0.0)] * pad)
    fh.write('float {}[] = {}\n'.format( '__' + name, '{'))
    fh.write('\t' + lines)
    fh.write('};\n')
    fh.write('_Mat {} = {}\n\t.nr = {},\n\t.nrq = {},\n\t.nc = {},\n\t.stride = {},\n\t.data.f = {}\n{};\n'.format('_' + name, '{', X.shape[0], nrq, 1, nrq * 4, '__' + name, '}'))
    fh.write('const flappie_matrix {} = &{};\n\n'.format(name, '_' + name))


def print_gru(gru, name):
    iW = _cudnn_to_guppy_gru(gru.cudnn_gru.weight_ih_l0)
    sW = _cudnn_to_guppy_gru(gru.cudnn_gru.weight_hh_l0)
    b = _cudnn_to_guppy_gru(gru.cudnn_gru.bias_ih_l0)
    cformatM(sys.stdout, '{}iW'.format(name), iW)
    cformatM(sys.stdout, '{}sW'.format(name), sW)
    cformatV(sys.stdout, '{}b'.format(name), b.reshape(-1))


def print_lstm(layer, name):
    iW = layer.lstm.weight_ih_l0
    sW = layer.lstm.weight_hh_l0
    b = layer.lstm.bias_ih_l0
    cformatM(sys.stdout, '{}iW'.format(name), iW)
    cformatM(sys.stdout, '{}sW'.format(name), sW)
    cformatV(sys.stdout, '{}b'.format(name), b.reshape(-1))


def print_convolution(convwrapper, name, scale=False):
    filterW =  convwrapper.conv.weight
    if scale:
        #  Scaling factor for MAD
        filterW *= 1.4826
    nfilter, nf , winlen = filterW.shape
    nf2 = 4 * math.ceil(nf / 4.0)
    cformatM(sys.stdout, '{}W'.format(name), filterW.permute(0, 2, 1).reshape(-1, nf),
             nr=nf2 * winlen - nf2 + nf, nc=nfilter)
    cformatV(sys.stdout, '{}b'.format(name), convwrapper.conv.bias.reshape(-1))
    sys.stdout.write("#define {}stride  {}\n".format(name, convwrapper.stride))
    sys.stdout.write("""#define {}nfilter  {}
    #define {}winlen  {}
    """.format(name, nfilter, name, winlen))





if __name__ == '__main__':
    args = parser.parse_args()
    modelid = args.id + '_'

    network = helpers.load_model(args.model)

    if isinstance(network.sublayers[0], DeltaSample):
        sys.stderr.write('* Removing initial DeltaSample layer\n')
        network.sublayers = network.sublayers[1:]

    sys.stdout.write("""#pragma once
    #ifndef FLIPFLOP_{}MODEL_H
    #define FLIPFLOP_{}MODEL_H
    #include "../util.h"
    """.format(modelid.upper(), modelid.upper()))

    """ Convolution layers
    """
    conv1 = network.sublayers[0]
    print_convolution(conv1, 'conv1_rnnrf_flipflop5_{}'.format(modelid),
                      scale=args.scale)
    conv2 = network.sublayers[1]
    print_convolution(conv2, 'conv2_rnnrf_flipflop5_{}'.format(modelid),
                      scale=args.scale)
    conv3 = network.sublayers[2]
    print_convolution(conv3, 'conv3_rnnrf_flipflop5_{}'.format(modelid),
                      scale=args.scale)

    """  Backward LSTM (first layer)
    """
    lstm1 = network.sublayers[3].layer
    print_lstm(lstm1, 'lstmB1_rnnrf_flipflop5_{}'.format(modelid))

    """  Forward LSTM (second layer)
    """
    lstm2 = network.sublayers[4]
    print_lstm(lstm2, 'lstmF2_rnnrf_flipflop5_{}'.format(modelid))

    """ backward LSTM(third layer)
    """
    lstm3 = network.sublayers[5].layer
    print_lstm(lstm3, 'lstmB3_rnnrf_flipflop5_{}'.format(modelid))

    """  Forward LSTM (fourth layer)
    """
    lstm4 = network.sublayers[6]
    print_lstm(lstm4, 'lstmF4_rnnrf_flipflop5_{}'.format(modelid))

    """ backward LSTM(fifth layer)
    """
    lstm5 = network.sublayers[7].layer
    print_lstm(lstm5, 'lstmB5_rnnrf_flipflop5_{}'.format(modelid))
    """ Global norm layer
    """
    gnlayer = network.sublayers[8]
    nstate = gnlayer.linear.weight.shape[0]
    cformatM(sys.stdout, 'FF_rnnrf_flipflop5_{}W'.format(modelid), gnlayer.linear.weight)
    cformatV(sys.stdout, 'FF_rnnrf_flipflop5_{}b'.format(modelid), gnlayer.linear.bias)

    sys.stdout.write('#endif /* FLIPFLOP_{}MODEL_H */'.format(modelid.upper()))
