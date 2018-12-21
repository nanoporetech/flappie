#!/usr/bin/env python
import argparse
from enum import Enum
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


BASE = "ACGTZ"
colour_scheme = { 'default' :
                          { 'A' : 'green', 'C' : 'blue', 'G' : 'orange', 'T' : 'red',
                            'Z' : 'purple', 'N' : 'grey'},
                  'eccles' :  # ACGT as suggested by David Eccles, http://www.gringene.org/
                          { 'A' : '#006400', 'C' : '#0000ff', 'G' : '#ffd700', 'T' : '#ff6347',
                            'Z' : '#cc79a7', 'N' : '#808080'},
                  'friendly' : # http://jfly.iam.u-tokyo.ac.jp/color/#pallet
                          { 'A' : '#009e73', 'C' : '#0072b2', 'G' : '#f0e442', 'T' : '#d55e00',
                            'Z' : '#cc79a7', 'N' : '#808080'},
                  'traditional' : 
                          { 'A' : 'green', 'C' : 'blue', 'G' : 'gold', 'T' : 'red',
                            'Z' : 'purple', 'N' : 'grey'}}

parser = argparse.ArgumentParser()
parser.add_argument('--analysis', metavar='number', type=int, default=0,
                    help='Analysis number for fast5 files')
parser.add_argument('--colours', '--colors', default='default',
                    choices=list(colour_scheme.keys()),
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


class FileType(Enum):
    flappie_trace = 0
    single_read_fast5 = 1
    multi_read_fast5 = 2



if __name__ == '__main__':
    args = parser.parse_args()

    colours = colour_scheme[args.colours]


    with h5py.File(args.hdf5, 'r') as h5:
        if 'file_version' in h5.attrs:
            if 'Raw' in h5:
                file_type = FileType.single_read_fast5
            else:
                file_type = FileType.multi_read_fast5
        else:
            file_type = FileType.flappie_trace
  
        if file_type == FileType.flappie_trace:
            #  Flappie format
            reads = list(h5.keys())
        else:
            #  Guppy fast5
            if file_type == FileType.single_read_fast5:
                #  Single-read fast5
                reads = [args.hdf5]
            else:
                #  Multi-read fast5
                reads = list(h5.keys())

        for read in reads:
            if file_type == FileType.flappie_trace:
               #  Flappie format
               sig = h5[posixpath.join(read, 'signal')][()]
               trace = h5[posixpath.join(read, 'trace')][()] / 255.0
            else:
               #  Guppy
               if file_type == FileType.single_read_fast5:
                   #  Guppy single-read
                   readh5 = h5
                   readno = list(readh5[posixpath.join('Raw', 'Reads')].keys())[0]
                   sig = readh5[posixpath.join('Raw', 'Reads', readno, 'Signal')][()] / 255.0
               else:
                   #  Guppy multi-read
                   readh5 = h5[read]
                   sig = readh5[posixpath.join('Raw', 'Signal')][()] / 255.0

               trace = readh5[posixpath.join('Analyses', 'Basecall_1D_{:03d}'.format(args.analysis),
                                             'BaseCalled_template', 'Trace')][()]
               segpath = posixpath.join('Analyses', 'Segmentation_{:03d}'.format(args.analysis),
                                        'Summary', 'segmentation')
               sig_start = readh5[segpath].attrs['first_sample_template']
               sig_length = readh5[segpath].attrs['duration_template']
               sig = sig[sig_start : sig_start + sig_length]
               
               
               
            nbase = trace.shape[1] // 2
            assert nbase * 2 == trace.shape[1]
            assert nbase == 4 or nbase == 5, "Unsupported number of bases"
            if args.flipflops:
                trace[:,nbase:] *= -1
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
            for i in range(nbase):
                pp.fill_between(x2, trace[:, i], color=colours[BASE[i]], alpha=0.3)
                pp.fill_between(x2, trace[:, i + nbase], color=colours[BASE[i]], alpha=0.3)

                pp.plot(x2, trace[:, i], color=colours[BASE[i]])
                pp.plot(x2, trace[:, i + nbase], color=colours[BASE[i]], linestyle='dashed')

            pp.grid()
            pp.show()
