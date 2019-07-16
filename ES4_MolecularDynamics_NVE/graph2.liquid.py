import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(15, 21))
fig.suptitle('Kripton, Liquid Phase', y='0.92', size='20', weight='20')
gs = gridspec.GridSpec(4, 2, figure=fig)

f4 = np.loadtxt('liquid/output_temp_2.dat')
ax2 = fig.add_subplot(gs[0, 0])
ax2.plot(f4)
ax2.set_xlabel("Number of simulation steps")
ax2.set_ylabel("Temperature (K)")
ax2.grid()
plt.title("Temperature trend")

x1, f1, error1 = np.loadtxt('liquid/output_ave_temp_2.dat', usecols=(0,1,2), unpack='true')
ax3= fig.add_subplot(gs[0, 1])
ax3.errorbar(x1, f1, yerr=error1, label='Temperatura media su ogni blocco')
ax3.set_xlabel("Number of blocks")
ax3.set_ylabel("Temperature (K)")
ax3.grid()
plt.title("Mean Temperature")

f1 = np.loadtxt('liquid/output_ekin_2.dat')
f2 = np.loadtxt('liquid/output_epot_2.dat')
f3 = np.loadtxt('liquid/output_etot_2.dat')
ax = fig.add_subplot(gs[1, :])
ax.plot(f1, label='Energia Cinetica')
ax.plot(f2, label='Energia Potenziale')
ax.plot(f3, label='Energia Totale')
ax.set_xlabel("Number of simulation steps")
ax.set_ylabel("Energy")
ax.grid()
ax.legend(loc='upper right')
plt.title('Energy trend')

x1, f1, error1 = np.loadtxt('liquid/output_ave_ekin_1.dat', usecols=(0,1,2), unpack='true')
x2, f2, error2 = np.loadtxt('liquid/output_ave_epot_1.dat', usecols=(0,1,2), unpack='true')
x3, f3, error3 = np.loadtxt('liquid/output_ave_etot_1.dat', usecols=(0,1,2), unpack='true')
ax6 = fig.add_subplot(gs[2, :])
ax6.errorbar(x1, f1, yerr=error1, label='Energia cinetica media su ogni blocco')
ax6.errorbar(x2, f2, yerr=error2, label='Energia potenziale media su ogni blocco')
ax6.errorbar(x3, f3, yerr=error3, label='Energia totale media su ogni blocco')
ax6.set_xlabel("Number of blocks")
ax6.set_ylabel("Energy")
ax6.grid()
ax6.legend(loc='right')
plt.title('Mean energies')

f4 = np.loadtxt('liquid/output_pres_2.dat')
ax4= fig.add_subplot(gs[3, 0])
ax4.plot(f4)
ax4.set_xlabel("Number of simulation steps")
ax4.set_ylabel("Pressure (Pa)")
ax4.grid()
plt.title("Pressure trend")

x1, f1, error1 = np.loadtxt('liquid/output_ave_pres_2.dat', usecols=(0,1,2), unpack='true')
ax5= fig.add_subplot(gs[3, 1])
ax5.errorbar(x1, f1, yerr=error1, label='Pressione media su ogni blocco')
#ax5.plot(x1, 0.8 +x1-x1)
ax5.set_xlabel("Number of blocks")
ax5.set_ylabel("Pressure (Pa)")
ax5.grid()
plt.title("Mean Pressure")

plt.show()
