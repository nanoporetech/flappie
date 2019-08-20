#!/usr/bin/env python3
import argparse
from itertools import islice
from multiprocessing import Pool
import numpy as np
import sys



class Positive(object):
    """Create an argparse argument type that accepts only positive values

    :param mytype: Type function for type to accept, e.g. `int` or `float`
    """

    def __init__(self, mytype):
        self.mytype = mytype

    def __repr__(self):
        return "positive {}".format(self.mytype)

    def __call__(self, y):
        yt = self.mytype(y)
        if yt <= 0:
            raise argparse.ArgumentTypeError('Argument must be {}'.format(self))
        return yt


parser = argparse.ArgumentParser()
parser.add_argument('--limit', default=None, type=Positive(int),
                    help='Limit number of reads processed')
parser.add_argument('--rlc', default=False, action='store_true',
                    help='Call run-length compressed sequence')
parser.add_argument('--no-rlc', dest='rlc', action='store_false',
                    help="Don't call run-length compressed sequence")
parser.add_argument('--run_max', default=50, type=Positive(int),
                    help='Maximum run for mean approximation')
parser.add_argument('--scale', default=(1.00, 1.00, 1.00, 1.00), nargs=4,
                    metavar=('scaleA', 'scaleC', 'scaleG', 'scaleT'), type=Positive(float),
                    help='Factors for per-base scale parameter')
parser.add_argument('--shape', default=(1.00, 1.00, 1.00, 1.00), nargs=4,
                    metavar=('shapeA', 'shapeC', 'shapeG', 'shapeT'), type=Positive(float),
                    help='Factors for per-base shape parameter')
parser.add_argument('--threads', default=1, type=Positive(int),
                    help='Number of threads to use')
parser.add_argument('--width', default=60, type=Positive(int),
                    help='Line width for Fasta output')
parser.add_argument('file', default='/dev/stdin', nargs='?')


def pow1p(x, y):
    return np.exp(y * np.log1p(x))


def run_estimate_mode(shape, scale, imax=50):
    inv_shape = np.reciprocal(shape)
    if shape <= 1.0:
        return 1
    run_mode = 1 + np.floor(scale * pow1p(-inv_shape, inv_shape))
    return run_mode.astype(int)


baseidx = {b : i for i, b in enumerate('ACGT')}
def runlength_basecall(read_data, shapef, scalef, imax=50):
    return ''.join([b * run_estimate_mode(sh * shapef[baseidx[b]],
                                          sc * scalef[baseidx[b]], imax=imax)
                    for b, sh, sc in read_data])


def read_generator(fh):
    first_read = True
    for line in fh:
        if line.startswith('#'):
            if not first_read:
                yield read_name, read_data
            #  New record
            first_read = False
            read_name = line[2:-1]
            read_data = []
        else:
            base, shape, scale = line.split('\t')
            read_data.append((base, float.fromhex(shape), float.fromhex(scale)))
    yield read_name, read_data


gbl_args = None
def init_basecall_worker(*args):
    global gbl_args
    if len(args) > 0:
        gbl_args = {'shape' : args[0], 'scale' : args[1], 'run_max' : args[2]}


def basecall_worker(indata):
    read_name, read_data = indata
    if gbl_args is None:
        basecall = ''.join([elt[0] for elt in read_data])
    else:
        basecall = runlength_basecall(read_data, gbl_args['shape'],
                                      gbl_args['scale'], gbl_args['run_max'])

    return read_name, basecall


if __name__ == '__main__':
    args = parser.parse_args()

    if args.rlc:
        init_params = []
    else:
        init_params = [args.shape, args.scale, args.run_max]

    with open(args.file, 'r') as fh:
        with Pool(processes=args.threads, initializer=init_basecall_worker, initargs=init_params) as pool:
            for res in pool.imap(basecall_worker, islice(read_generator(fh), args.limit)):
                read_name, basecall = res
                print('>{}'.format(read_name))
                print('\n'.join([basecall[st : st + args.width]
                                 for st in range(0, len(basecall), args.width)]))
