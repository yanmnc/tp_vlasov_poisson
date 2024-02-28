import numpy as np

def fourier_2d(F0, y0, x0):
    nx0 = len(x0)
    nx = 2 * (nx0 // 2)
    hnx = nx // 2
    ny0 = len(y0)
    ny = 2 * (ny0 // 2)
    hny = ny // 2

    x = x0[:nx]
    y = y0[:ny]
    F = F0[:ny, :nx]

    # Creation of the vector in the Fourier space ky(1:ny) kx(1:nx)
    Lx = x[-1] - x[0]
    dx = x[1] - x[0]
    dkx = 2 * np.pi / (Lx + dx)
    kx = np.zeros(nx)
    kx[:hnx] = -(np.arange(hnx, 0, -1)) * dkx
    kx[hnx:nx] = np.arange(0, hnx) * dkx

    Ly = y[-1] - y[0]
    dy = y[1] - y[0]
    dky = 2 * np.pi / (Ly + dy)
    ky = np.zeros(ny)
    ky[:hny] = -(np.arange(hny, 0, -1)) * dky
    ky[hny:ny] = np.arange(0, hny) * dky

    # Fourier transformation using numpy's fft2
    TFF = np.zeros_like(F, dtype=np.complex128)
    var = np.fft.fft2(F) / (nx * ny)
    AA = np.zeros_like(var)
    AA[:, :hnx] = var[:, hnx:nx]
    AA[:, hnx:nx] = var[:, :hnx]
    TFF[:hny, :] = AA[hny:ny, :]
    TFF[hny:ny, :] = AA[:hny, :]

    return TFF, ky, kx

# Example usage:
# TFF, ky, kx = fourier_2d(F0, y0, x0)
