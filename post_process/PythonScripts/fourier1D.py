import numpy as np

def fourier_1d(F0, x0):
    nx0 = len(x0)
    nx = 2 * (nx0 // 2)
    hnx = nx // 2
    x = x0[:nx]
    F = F0[:nx]

    # Creation of the vector in the Fourier space kx(1:nx)
    Lx = x[-1] - x[0]
    dx = x[1] - x[0]
    dkx = 2 * np.pi / (Lx + dx)
    kx = np.zeros(nx)
    kx[:hnx] = -(np.arange(hnx, 0, -1)) * dkx
    kx[hnx:nx] = np.arange(0, hnx) * dkx

    # Fourier transformation using numpy's fft
    TFF = np.zeros_like(F, dtype=np.complex128)
    var = np.fft.fft(F) / nx
    TFF[:hnx] = var[hnx:nx]
    TFF[hnx:nx] = var[:hnx]

    return TFF, kx