import numpy as np
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, r2_score
import argparse

# This script calculates R2 and other metrics

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-est", "--est", help = "Path to y estimates file")
parser.add_argument("-true", "--true", help = "Path to true phen file")
parser.add_argument("-model", "--model", help = "Determine wheter probit or linear")
args = parser.parse_args()
est_fpath = args.est
phen_fpath = args.true
model = args.model

print(est_fpath)
print(phen_fpath)

def load_file(path, col):
    y = []
    file = open(path, "rb")
    for row in file:
        row = row.split()
        y.append(float(row[col]))
    return np.array(y)

y_est = load_file(est_fpath, 0)
y_true = load_file(phen_fpath, 2)

r2 = r2_score(y_true, y_est)

print("R2 = %0.4f" % ( r2))

if model == "probit":
    y_est[y_est > 0] = 1
    y_est[y_est < 0] = 0

    acc = accuracy_score(y_true, y_est)
    fscore = f1_score(y_true, y_est)
    tn, fp, fn, tp = confusion_matrix(y_true, y_est).ravel()
    fdr = fp / (fp + tp)

    print("ACC = %0.4f, F score = %0.4f, FDR = %0.4f, R2 = %0.4f" % (acc, fscore, fdr, r2))
    print("TP: %d, TN: %d, FP: %d, FN: %d" % (tp, tn, fp, fn))