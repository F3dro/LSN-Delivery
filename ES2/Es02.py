import numpy as np
import matplotlib
import matplotlib.pyplot as plt

x, f, error = np.loadtxt('Es02.1.1.res', usecols=(0,1,2), unpack='true')
plt.errorbar(x, f, yerr=error)
plt.xlabel("Number of blocks")
plt.ylabel("f(x)")
plt.grid()
plt.show()
