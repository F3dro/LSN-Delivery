import numpy as np
import matplotlib
import matplotlib.pyplot as plt

x=np.arange(0.0, 1.0, 0.015)
f1=np.pi/2*np.cos(np.pi/2*x)
f2=-2*x+2
f3= x-x+1
fig, ax = plt.subplots()
ax.plot(x, f1, label=r'$f(x)=\frac{\pi}{2}\cdot cos(\frac{\pi}{2}\cdot x)$')
ax.plot(x, f2, label=r'$f(x)=2-2x$')
ax.plot(x, f3, label=r'$y=1$')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Scelta della funzione per l\'importance sampling')
plt.grid()
legend=ax.legend()
plt.show()
