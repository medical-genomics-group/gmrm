import matplotlib.pyplot as plt
import numpy as np
import csv
import argparse

# This script plot the MCMC traces for noise and genetic variance, number of included markers, and heritability

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-path", "--path", help = "Path to csv file")
parser.add_argument("-iterations", "--iterations", help = "Iterations range start:end")
args = parser.parse_args()
path = args.path
ran = args.iterations
start = int(ran.split(':')[0])
end = int(ran.split(':')[1])

sigmag = []
sigmae = []
h2 = []
mincl = []
with open(path) as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        #print(row)
        sigmag.append(float(row[2]))
        sigmae.append(float(row[3]))
        h2.append(float(row[4]))
        mincl.append(float(row[5]))

sigmag = np.array(sigmag[start:end])
sigmae = np.array(sigmae[start:end])
h2 = np.array(h2[start:end])
mincl = np.array(mincl[start:end])

print("Means: ")
print("h2 = %0.4f" % np.mean(h2))
print("Incl. markers = %d" % np.mean(mincl))

fig, ax = plt.subplots(4,2)
fig.tight_layout(pad=1)

ax[0,0].set_title("Genetic variance $\u03C3_G$")
ax[0,0].plot(sigmag)
ax[0,0].set_xlabel("Iterations")
ax[0,0].set_ylabel("$\u03C3_G$")

ax[0,1].hist(sigmag, bins=20)
ax[0,1].set_xlabel("$\u03C3_G$")
ax[0,1].set_ylabel("Frequency")

ax[1,0].set_title("Noise variance $\u03C3_E$")
ax[1,0].plot(sigmae)

ax[1,0].set_xlabel("Iterations")
ax[1,0].set_ylabel("$\u03C3_E$")

ax[1,1].hist(sigmae, bins=20)
ax[1,1].set_xlabel("$\u03C3_E$")
ax[1,1].set_ylabel("Frequency")

ax[3,0].set_title("Number of included markers")
ax[3,0].plot(mincl)
ax[3,0].set_xlabel("Iterations")
ax[3,0].set_ylabel("Incl. m.")

ax[3,1].hist(mincl, bins=20)
ax[3,1].set_xlabel("Included markers")
ax[3,1].set_ylabel("Frequency")

ax[2,0].set_title("Heritability $h^2$")
ax[2,0].plot(h2)
ax[2,0].set_xlabel("Iterations")
ax[2,0].set_ylabel("$h^2$")

ax[2,1].hist(h2, bins=20, density=True)
ax[2,1].set_xlabel("$h^2$")
ax[2,1].set_ylabel("Frequency")


plt.show()