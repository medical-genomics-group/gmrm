import numpy as np
import argparse
import struct
from array import array
import os

# This script calculate posterior inclusion probabilities

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-bet", "--bet", help = "Path to beta bin file")
parser.add_argument("-iterations", "--iterations", help = "Number of iterations")
parser.add_argument("-M", "--M", help = "Number of markers")
args = parser.parse_args()
betfile = args.bet
iterations = args.iterations
M = int(args.M)
it_start = int(iterations.split(":")[0])
it_end = int(iterations.split(":")[1])
iterations = (it_end - it_start)

basename = os.path.basename(betfile)
basename = basename.split('.')[0]
dirpath = os.path.dirname(betfile)
print(betfile)
print(basename)
print(dirpath)

pip = np.zeros(M)
with open(betfile, "rb") as f:
    buffer = f.read(4)
    [m] = struct.unpack('I', buffer)
    print("Number of markers: ", m)
    for i in range(it_end):
        buffer = f.read(4)
        [it] = struct.unpack('I', buffer)
        buffer = f.read(m*8)
        if it >= it_start:
            beta = struct.unpack(str(m)+'d', buffer)
            beta = np.array(beta)
            beta[np.abs(beta) > 0] = 1
            pip += beta 
pip /= iterations

output_file = open(os.path.join(dirpath, basename+'.pip'), 'wb')
float_array = array('d', pip)
float_array.tofile(output_file)
output_file.close()