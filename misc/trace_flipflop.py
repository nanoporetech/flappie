#!/usr/bin/env python
import argparse
import h5py
from matplotlib import pyplot as pp
import numpy as np
import posixpath


class AutoBool(argparse.Action):

    def __init__(self, option_strings, dest, default=None, required=False, help=None):
        """Automagically create --foo / --no-foo argument pairs"""

        if default is None:
            raise ValueError('You must provide a default with AutoBool action')
        if len(option_strings) != 1:
            raise ValueError('Only single argument is allowed with AutoBool action')
        opt = option_strings[0]
        if not opt.startswith('--'):
            raise ValueError('AutoBool arguments must be prefixed with --')

        opt = opt[2:]
        opts = ['--' + opt, '--no-' + opt]
        if default:
            default_opt = opts[0]
        else:
            default_opt = opts[1]
        super(AutoBool, self).__init__(opts, dest, nargs=0, const=None,
                                       default=default, required=required,
                                       help='{} (Default: {})'.format(help, default_opt))

    def __call__(self, parser, namespace, values, option_strings=None):
        if option_strings.startswith('--no-'):
            setattr(namespace, self.dest, False)
        else:
            setattr(namespace, self.dest, True)

    @staticmethod
    def filter_option_strings(strings):
        for s in strings:
            s = s.strip('-')
            if s[:3] != 'no-':
                yield s


class FileExists(argparse.Action):
    """Check if the input file exist."""

    def __call__(self, parser, namespace, values, option_string=None):
        if not os.path.exists(values):
            raise RuntimeError("File/path for '{}' does not exist, {}".format(self.dest, values))
        setattr(namespace, self.dest, values)


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


class Maybe(object):
    """Create an argparse argument type that accepts either given type or 'None'

    :param mytype: Type function for type to accept, e.g. `int` or `float`
    """

    def __init__(self, mytype):
        self.mytype = mytype

    def __repr__(self):
        return "None or {}".format(self.mytype)

    def __call__(self, y):
        try:
            if y == 'None':
                res = None
            else:
                res = self.mytype(y)
        except:
            raise argparse.ArgumentTypeError('Argument must be {}'.format(self))
        return res


BASE = "ACGT"
colour_scheme = { 'default' :
                          { 'A' : 'green', 'C' : 'blue', 'G' : 'orange', 'T' : 'red'},
                  'friendly' : # http://jfly.iam.u-tokyo.ac.jp/color/#pallet
                          { 'A' : '#003c32', 'C' : '#002d46', 'G' : 'grey', 'T' : '#502800'},
                  'gringer' :
                          { 'A' : '#006400', 'C' : '#0000ff', 'G' : '#ffd700', 'T' : '#ff6347'},
                  'traditional' : 
                          { 'A' : 'green', 'C' : 'blue', 'G' : 'grey', 'T' : 'red'}}

parser = argparse.ArgumentParser()
parser.add_argument('--colours', default='default', choices=list(colour_scheme.keys()),
                    help='Change trace colour scheme')
parser.add_argument('--depop', default=None, metavar='threshold',
                    type=Maybe(Positive(float)), help='Filter pops from signal')
parser.add_argument('--limit', default=10, type=Maybe(int))
parser.add_argument('--flipflops', default=False, action=AutoBool,
                    help='Plot the flop states as negative probabilites')
parser.add_argument('hdf5')



def depop(sig, thresh):
    where_big = np.where(np.abs(sig) > thresh)[0]
    sig[where_big] = 0
    return sig



if __name__ == '__main__':
    args = parser.parse_args()

    colours = colour_scheme[args.colours]


    with h5py.File(args.hdf5, 'r') as h5:
        reads = list(h5.keys())

        for read in reads:
            sig = h5[posixpath.join(read, 'signal')][()]
            trace = h5[posixpath.join(read, 'trace')][()] / 255.0
            assert trace.shape[1] == 8, "Trace viewer does not yet support modified bases"
            if args.flipflops:
                trace[:,4:] *= -1
            down_sample_factor = round(len(sig) / float(len(trace)))

            if args.depop is not None:
                sig = depop(sig, args.depop)

            ax1 = pp.subplot(211)
            pp.title(read)
            pp.ylabel('Normalised signal')
            pp.plot(np.arange(len(sig)), sig, color='grey')

            pp.subplot(212, sharex=ax1)
            pp.xlabel('time (samples)')
            pp.ylabel('State probability')

            x2 = down_sample_factor * np.arange(len(trace))
            pp.fill_between(x2, trace[:,0], color=colours['A'], alpha=0.3)
            pp.fill_between(x2, trace[:,1], color=colours['C'], alpha=0.3)
            pp.fill_between(x2, trace[:,2], color=colours['G'], alpha=0.3)
            pp.fill_between(x2, trace[:,3], color=colours['T'], alpha=0.3)
            pp.fill_between(x2, trace[:,4], color=colours['A'], alpha=0.3)
            pp.fill_between(x2, trace[:,5], color=colours['C'], alpha=0.3)
            pp.fill_between(x2, trace[:,6], color=colours['G'], alpha=0.3)
            pp.fill_between(x2, trace[:,7], color=colours['T'], alpha=0.3)

            pp.plot(x2, trace[:,0], color=colours['A'])
            pp.plot(x2, trace[:,1], color=colours['C'])
            pp.plot(x2, trace[:,2], color=colours['G'])
            pp.plot(x2, trace[:,3], color=colours['T'])
            pp.plot(x2, trace[:,4], color=colours['A'], linestyle='dashed')
            pp.plot(x2, trace[:,5], color=colours['C'], linestyle='dashed')
            pp.plot(x2, trace[:,6], color=colours['G'], linestyle='dashed')
            pp.plot(x2, trace[:,7], color=colours['T'], linestyle='dashed')

            pp.grid()
            pp.show()
