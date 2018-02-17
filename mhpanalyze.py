import numpy as np
import sys

filename = sys.argv[1]
ids = [int(s) for s in sys.argv[2].split(',')]
beta_vals = np.zeros((1, len(ids)+1))
row = 0
col = 1
with open(filename) as fp:
    for i, line in enumerate(fp):
        check = [i for i in ['REMARK', 'END'] if i in line]
        if not check:
            index = int(line.split()[1])
            if index in ids:
                print(row, col)
                beta_vals[row][col] = float(line.split()[9])
                print(beta_vals)
                col += 1
                if col == len(ids) + 1:
                    col = 1
                    row += 1
                    beta_vals = np.vstack((beta_vals, np.zeros((1, len(ids)+1))))
                    beta_vals[row][0] = row

print(beta_vals)
