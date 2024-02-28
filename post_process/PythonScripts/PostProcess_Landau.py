import numpy as np
import matplotlib.pyplot as plt

from read_VlasovPoiss import read_VlasovPoiss
from fourier1D import fourier_1d

s = read_VlasovPoiss()
kx0 = s['kx0']
Enpot = s['Enpot']
Enkin = s['Enkin']
time = s['time']
entropy = s['entropy']
Phi1D_evol = s['Phi1D_evol']
Ntime = s['Ntime']
vg = s['vg']
f2D_evol_true = s['f2D_evol_true']
f1Dv_evol = s['f1Dv_evol']

# Approximate solution for k<<1
omega_th = np.sqrt(1 + np.sqrt(1 + 12 * kx0**2)) / np.sqrt(2)
gamma_th = -np.sqrt(np.pi / 2) * omega_th / abs(kx0) * np.exp(-omega_th**2 / (2 * kx0**2)) / (2 * kx0**2 / omega_th**3 * (1 + 6 * kx0**2 / omega_th**2))
vphi_th = omega_th / kx0

print('*********************************')
print('Wave vector k=', kx0)
print('Using "DispRelation_Landau.m" to solve the dispersion relation, enter:')
omega_str = input('   Real frequency [default: limit k<<1]:      omega = ')
gamma_str = input('   Imaginary frequency [default: limit k<<1]: gamma = ')
print('*********************************')
if omega_str:
    omega_th = float(omega_str)
if gamma_str:
    gamma_th = float(gamma_str)
vphi_th = np.sqrt(omega_th**2 + gamma_th**2) / kx0

# Plot the potential energy
Enpot_th = Enpot[0] * np.cos(omega_th * time)**2 * np.exp(2 * gamma_th * time)

plt.figure(1)
plt.semilogy(time, Enpot, 'r+-', time, Enpot_th, 'b--')
plt.grid(True)
plt.xlabel('time')
plt.ylabel('Potential energy')
plt.legend(['Enpot', 'Enpot_th'])
plt.show()

# Plot energy & entropy conservation
Entot = Enkin + Enpot
Entot_init = (Enkin[0] + Enpot[0]) / 2.

plt.figure(2)
plt.subplot(211)
plt.plot(time, Enkin, 'r+-', time, Enpot, 'b-o')
plt.plot(time, Entot-Entot_init, 'g+-')
plt.grid(True)
plt.xlabel('time')
plt.ylabel('Energy')
plt.legend(['Enkin', 'Enpot', 'Entot-Entot_init'])

plt.subplot(212)
plt.semilogy(time, (entropy - entropy[0]) / entropy, 'r+-')
plt.grid(True)
plt.xlabel('time')
plt.ylabel('Entropy')
plt.show()

# Real frequency - comparison to theoretical value
FTPhi_t, omega = fourier_1d(Phi1D_evol[0, :], time)
iomega_max = np.argmax(np.abs(FTPhi_t))
omega_max = omega[iomega_max]

print('=========================')
print('    omega at max[FTPhi] = ', omega_max)
print('=========================')

plt.figure(3)
plt.subplot(211)
plt.plot(time, Phi1D_evol[0, :], '-r.')
plt.grid(True)
plt.xlabel('time')
plt.ylabel('Phi(t,x=0)')

plt.subplot(212)
plt.semilogy(omega, np.abs(FTPhi_t), '-r.')
plt.grid(True)
plt.xlabel('frequency')
plt.ylabel('|FT[phi(t,x=0)]|')
plt.plot(omega_max * np.ones(2), [min(np.abs(FTPhi_t)), max(np.abs(FTPhi_t))], 'r--')
plt.plot(omega_th * np.ones(2), [min(np.abs(FTPhi_t)), max(np.abs(FTPhi_t))], 'k--', linewidth=2)
plt.legend(['Spectrum', 'max of spectrum', 'theoretical frequency'])
plt.show()

# Small scales in velocity space
dfv_evol = f1Dv_evol - np.outer(f1Dv_evol[:, 0], np.ones(Ntime))

FTdfv = np.zeros((Ntime, len(vg)), dtype=complex)
for it in range(Ntime):
    FTdfv[it, :], kv = fourier_1d(dfv_evol[:, it], vg)

dkv = kv[1] - kv[0]

plt.figure(4)
plt.pcolormesh(kv - dkv / 2, time, np.abs(FTdfv), shading='auto')
plt.colorbar()
plt.xlabel('k_v')
plt.ylabel('time')
plt.plot(kx0 * time + 2 * dkv, time, 'k--')
plt.plot(-kx0 * time - 2 * dkv, time, 'k--')
plt.title('Fourier Transform of delta f in velocity space')
plt.show()

# f2D plot
if f2D_evol_true:
    Nx = s['Nx']
    xg = s['xg']
    f2D_evol = s['f2D_evol']
    fM = np.exp(-vg**2 / 2) / np.sqrt(2 * np.pi)
    fM2D = np.outer(np.ones(Nx), fM)

    for itdiag in range(0, Ntime, 5):
        plt.figure(5)
        plt.pcolormesh(xg, vg, (f2D_evol[:, :, itdiag] - f2D_evol[:, :, 0]), shading='auto')
        plt.colorbar()
        plt.plot([xg[0], xg[-1]], [vphi_th, vphi_th], 'r', label='vphi_th')
        plt.xlabel('x coordinate')
        plt.ylabel('velocity')
        plt.title('delta f(x,v) at time ' + str(time[itdiag]) + ' / ' + str(time[-1]))
        plt.legend()
        plt.pause(0.1)
        plt.show()
