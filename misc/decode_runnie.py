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
parser.add_argument('--run_max', default=50, type=Positive(int),
                    help='Maximum run for mean approximation')
parser.add_argument('--scaling_factor', default=(1.0, 1.05), nargs=2,
                    metavar=('shape', 'scale'), type=Positive(float),
                    help='Scaling for parameters')
parser.add_argument('--threads', default=1, type=Positive(int),
                    help='Number of threads to use')
parser.add_argument('--width', default=60, type=Positive(int),
                    help='Line width for Fasta output')
parser.add_argument('file', default='/dev/stdin', nargs='?')


def cp_weibull(x, shape, scale):
    """  Tail cumulative probability of Weibull distribution

        :param x: points at which to evaluate
        :param shape: Shape parameter
        :param scale: Scale parameter
    """
    return np.exp(-np.power(x / scale, shape))


def discrete_weibull(i, shape, scale):
    """  Probability mass for Discrete Weibull distribution

        :param i: points at which to evaluate
        :param shape: Shape parameter
        :param scale: Scale parameter
    """
    return cp_weibull(i, shape, scale) - cp_weibull(i + 1, shape, scale)


def approx_mean_discrete_weibull(imax, shape, scale):
    """  Approximate mean of Discrete Weibull

        Uses E(X) = \sum P(X>x)

        :param imax: Upper limit to evaluate at
        :param shape: Shape parameter
        :param scale: Scale parameter
    """
    i = np.arange(1, imax, dtype='f4')
    return np.sum(cp_weibull(i, shape, scale))


def run_estimate(shape, scale, imax=50):
    return 1 + np.round(approx_mean_discrete_weibull(imax, shape, scale)).astype(int)


def runlength_basecall(read_data, shapef, scalef, imax=50):
    return ''.join([b * run_estimate(sh * shapef, sc * scalef, imax=imax) for b, sh, sc in read_data])


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


gbl_shape_factor = None
gbl_scale_factor = None
gbl_run_max = None
def init_basecall_worker(shape_factor, scale_factor, run_max):
    global gbl_shape_factor, gbl_scale_factor, gbl_run_max
    gbl_shape_factor = shape_factor
    gbl_scale_factor = scale_factor
    gbl_run_max = run_max


def basecall_worker(indata):
    read_name, read_data = indata
    basecall = runlength_basecall(read_data, gbl_shape_factor, gbl_scale_factor, gbl_run_max)
    
    return read_name, basecall


if __name__ == '__main__':
    args = parser.parse_args()

    init_params = [args.scaling_factor[0], args.scaling_factor[1], args.run_max]

    with open(args.file, 'r') as fh:
        with Pool(processes=args.threads, initializer=init_basecall_worker, initargs=init_params) as pool:
            for res in pool.imap(basecall_worker, islice(read_generator(fh), args.limit)):
                read_name, basecall = res
                print('>{}'.format(read_name))
                print('\n'.join([basecall[st : st + args.width]
                                 for st in range(0, len(basecall), args.width)]))
