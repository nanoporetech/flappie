#  Copyright 2018 Oxford Nanopore Technologies, Ltd

#  This Source Code Form is subject to the terms of the Oxford Nanopore
#  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
#  was not distributed with this file, You can obtain one at
#  http://nanoporetech.com


import argparse
import numpy as np

parser = argparse.ArgumentParser(description='')
parser.add_argument('file')

def read_crp(filename):
    with open(filename, 'r') as fh:
        nr, nc = [int(x) for x in fh.readline().rstrip().split()]

        mat = np.zeros((nc, nr))
        for col in range(nc):
            mat[col] = [float.fromhex(x) for x in fh.readline().rstrip().split()]

    return mat

if __name__ == '__main__':
    args = parser.parse_args()

    mat = read_crp(args.file)
    print(mat)
