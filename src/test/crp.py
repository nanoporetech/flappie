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
