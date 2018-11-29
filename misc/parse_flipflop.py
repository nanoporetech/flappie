#!/usr/bin/env python3

#  Copyright 2018 Oxford Nanopore Technologies, Ltd

#  This Source Code Form is subject to the terms of the Oxford Nanopore
#  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
#  was not  distributed with this file, You can obtain one at
#  http://nanoporetech.com

import argparse
import pickle
import math
import numpy as np
import re
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--id', default='' , help='Identifier for model names')
parser.add_argument('model', help='Pickle to read model from')


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


if __name__ == '__main__':
    args = parser.parse_args()
    modelid = args.id + '_'

    with open(args.model, 'rb') as fh:
        network = pickle.load(fh, encoding='latin1')
    network_major_version = network.version[0] if isinstance(network.version, tuple) else network.version
    assert network_major_version >= 2, "Sloika model must be version >= 2 but model is {}.\nPerhaps you need to run Sloika's model_upgrade.py".format(network.version)

    sys.stdout.write("""#pragma once
    #ifndef FLIPFLOP_{}MODEL_H
    #define FLIPFLOP_{}MODEL_H
    #include "../util.h"
    """.format(modelid.upper(), modelid.upper()))

    """ Convolution layer
    """

    filterW =  network.sublayers[0].W.get_value()
    nfilter, _ , winlen = filterW.shape
    cformatM(sys.stdout, 'conv_rnnrf_flipflop_{}W'.format(modelid), filterW.reshape(-1, 1), nr = winlen * 4 - 3, nc=nfilter)
    cformatV(sys.stdout, 'conv_rnnrf_flipflop_{}b'.format(modelid), network.sublayers[0].b.get_value().reshape(-1))
    sys.stdout.write("#define conv_rnnrf_flipflop_{}stride  {}\n".format(modelid, network.sublayers[0].stride))
    sys.stdout.write("""#define {}nfilter  {}
    #define _conv_rnnrf_flipflop_{}winlen  {}
    """.format(modelid, nfilter, modelid, winlen))

    """  Backward GRU (first layer)
    """
    gru1 = network.sublayers[1].sublayers[0].sublayers[0]
    cformatM(sys.stdout, 'gruB1_rnnrf_flipflop_{}iW'.format(modelid), gru1.iW.get_value())
    cformatM(sys.stdout, 'gruB1_rnnrf_flipflop_{}sW'.format(modelid), gru1.sW.get_value())
    cformatM(sys.stdout, 'gruB1_rnnrf_flipflop_{}sW2'.format(modelid), gru1.sW2.get_value())
    cformatV(sys.stdout, 'gruB1_rnnrf_flipflop_{}b'.format(modelid), gru1.b.get_value().reshape(-1))

    """  Forward GRU (second layer)
    """
    gru2 = network.sublayers[2].sublayers[0]
    cformatM(sys.stdout, 'gruF2_rnnrf_flipflop_{}iW'.format(modelid), gru2.iW.get_value())
    cformatM(sys.stdout, 'gruF2_rnnrf_flipflop_{}sW'.format(modelid), gru2.sW.get_value())
    cformatM(sys.stdout, 'gruF2_rnnrf_flipflop_{}sW2'.format(modelid), gru2.sW2.get_value())
    cformatV(sys.stdout, 'gruF2_rnnrf_flipflop_{}b'.format(modelid), gru2.b.get_value().reshape(-1))

    """ backward GRU(third layer)
    """
    gru3 = network.sublayers[3].sublayers[0].sublayers[0]
    cformatM(sys.stdout, 'gruB3_rnnrf_flipflop_{}iW'.format(modelid), gru3.iW.get_value())
    cformatM(sys.stdout, 'gruB3_rnnrf_flipflop_{}sW'.format(modelid), gru3.sW.get_value())
    cformatM(sys.stdout, 'gruB3_rnnrf_flipflop_{}sW2'.format(modelid), gru3.sW2.get_value())
    cformatV(sys.stdout, 'gruB3_rnnrf_flipflop_{}b'.format(modelid), gru3.b.get_value().reshape(-1))

    """  Forward GRU (fourth layer)
    """
    gru4 = network.sublayers[4].sublayers[0]
    cformatM(sys.stdout, 'gruF4_rnnrf_flipflop_{}iW'.format(modelid), gru4.iW.get_value())
    cformatM(sys.stdout, 'gruF4_rnnrf_flipflop_{}sW'.format(modelid), gru4.sW.get_value())
    cformatM(sys.stdout, 'gruF4_rnnrf_flipflop_{}sW2'.format(modelid), gru4.sW2.get_value())
    cformatV(sys.stdout, 'gruF4_rnnrf_flipflop_{}b'.format(modelid), gru4.b.get_value().reshape(-1))

    """ backward GRU(fifth layer)
    """
    gru5 = network.sublayers[5].sublayers[0].sublayers[0]
    cformatM(sys.stdout, 'gruB5_rnnrf_flipflop_{}iW'.format(modelid), gru5.iW.get_value())
    cformatM(sys.stdout, 'gruB5_rnnrf_flipflop_{}sW'.format(modelid), gru5.sW.get_value())
    cformatM(sys.stdout, 'gruB5_rnnrf_flipflop_{}sW2'.format(modelid), gru5.sW2.get_value())
    cformatV(sys.stdout, 'gruB5_rnnrf_flipflop_{}b'.format(modelid), gru5.b.get_value().reshape(-1))
    """ Global norm layer
    """
    nstate = network.sublayers[6].W.get_value().shape[0]
    cformatM(sys.stdout, 'FF_rnnrf_flipflop_{}W'.format(modelid), network.sublayers[6].W.get_value())
    cformatV(sys.stdout, 'FF_rnnrf_flipflop_{}b'.format(modelid), network.sublayers[6].b.get_value())

    sys.stdout.write('#endif /* FLIPFLOP_{}MODEL_H */'.format(modelid.upper()))
