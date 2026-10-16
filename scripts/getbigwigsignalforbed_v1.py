import sys
cores = int(sys.argv[1])
import os
import pyBigWig
import numpy as np

import sys

if len(sys.argv) != 4:
    print("usage: python", sys.argv[0], "bw_file bed_file", file=sys.stderr)
    sys.exit(1)  

bw_file = sys.argv[2]
bed_file = sys.argv[3]

if not os.path.exists(bw_file):
    print(f"{bw_file} does not exist.", file=sys.stderr)
    sys.exit(1)  

if not os.path.exists(bed_file):
    print(f"{bed_file} does not exist.", file=sys.stderr)
    sys.exit(1)  

print(f"# Input bigWig file: {bw_file}")
print(f"# Input bed file: {bed_file}")

bw = pyBigWig.open(bw_file)

print("chr", "start", "end", "my_mean", sep="\t")

for line in open(bed_file):
    cols = line.strip().split()
    vals = bw.values(cols[0], int(cols[1]), int(cols[2]))
    vals1 = np.nan_to_num(vals, nan=0.0)
    my_mean = np.mean(vals1)
    print(cols[0], int(cols[1]), int(cols[2]), my_mean, sep="\t")
bw.close()
