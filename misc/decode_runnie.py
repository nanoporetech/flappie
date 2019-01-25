#!/usr/bin/env python
import argparse
from itertools import chain, islice
import numpy as np
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--width', default=60, type=int)
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
    i = np.arange(1, imax)
    return np.sum(cp_weibull(i, shape, scale))


def run_estimate(shape, scale, imax=10):
    return 1 + np.round(approx_mean_discrete_weibull(imax, shape, scale)).astype(int)


def runlength_basecall(read_data):
    return ''.join([b * run_estimate(sh, sc) for b, sh, sc in read_data])


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


if __name__ == '__main__':
    args = parser.parse_args()

    with open(args.file, 'r') as fh:
        for read_name, read_data in read_generator(fh):
            basecall = runlength_basecall(read_data)
            print('>{}'.format(read_name))
            print('\n'.join([basecall[st : st + args.width]
                             for st in range(0, len(basecall), args.width)]))
