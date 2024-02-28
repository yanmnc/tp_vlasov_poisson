import os
import numpy as np

from create_file_list import create_file_list


def read_1dfile_evol(prefix, name):
    '''creates a file list from the acronym/name files
    then aggregates the result and returns it'''

    print(f'--> Reading of the {name}_<num>.dat files')
    acronym = prefix + name
    filelist, nbfiles = create_file_list(acronym, '.dat', 0)
    tmp = np.loadtxt(filelist[0])

    # aggregates 1D fields over time
    evol = np.empty((len(tmp[:, 1]), nbfiles))  
    for ind in range(nbfiles):
        tmp = np.loadtxt(filelist[ind])
        evol[:, ind] = tmp[:, 1]
    
    lenfield = tmp.shape[0]
    lenaggregate = nbfiles
    grid = tmp[:, 0]
    
    return evol, lenfield, lenaggregate, grid


def read_VlasovPoiss():
    # Choice of the directory (if necessary)
    directory = input('Directory for data reading [default = current] ?: ')
    prefix = directory + '/' if directory else ''
    
    nb_s = 0  # Counter for the number of variables saved in the structure s

    # Read the file 'param_simu.dat'
    print('--> Reading of the param_simu.dat file')
    f_param = np.loadtxt(prefix + 'param_simu.dat')
    name = ['kx0', 'perturb', 'v0', 'T0', 'epsilon', 'vdiag_forf1Dx', 'nu_coll', 'Phi_ext0', 'kx_ext', 'omega_ext', 'coef_OhmsLaw', 'E0']
    nbparam = len(name)
    s = {}
    for i in range(nbparam):
        s[name[i]] = f_param[i]

    # Read the file 'VlasovPoiss_res.dat'
    print('--> Reading of the VlasovPoiss_res.dat file')
    f_pb2d = np.loadtxt(prefix + 'VlasovPoiss_res.dat')
    names = ['time', 'nbions', 'Enkin', 'Enpot', 'entropy', 'L1_norm', 'L2_norm']
    for ind, name in enumerate(names):
        s[name] = f_pb2d[:, ind]
    
    # Read all the files 'Phi1D_<num>.dat'
    Phi1D_evol, Nx, Ntime, xg = read_1dfile_evol(prefix, 'Phi1D')
    
    # Read all the files 'dens_<num>.dat'
    dens_evol, _, _, _ = read_1dfile_evol(prefix, 'dens')
    
    # Read all the files 'velo_<num>.dat'
    velo_evol, _, _, _ = read_1dfile_evol(prefix, 'velo')

    # Read all the files 'temp_<num>.dat'
    temp_evol, _, _, _ = read_1dfile_evol(prefix, 'temp')

    # Read all the files 'flux_<num>.dat'
    flux_evol, _, _, _ = read_1dfile_evol(prefix, 'flux')
    
    # Read all the files 'f1Dx_<num>.dat'
    f1Dx_evol, _, _, _ = read_1dfile_evol(prefix, 'f1Dx')

    # Read all the files 'f1Dv_<num>.dat'
    f1Dv_evol, Nv, _, vg = read_1dfile_evol(prefix, 'f1Dv')
    
    # Read all the files 'f2D_<num>.dat'
    f2D_acronym = prefix + 'f2D'
    f2D_filelist, f2D_nbfiles = create_file_list(f2D_acronym, '.dat', 0)

    if f2D_nbfiles > 1:
        print(f'--> Reading of the f2D_<num>.dat files')
        f2D_evol = np.empty((Nv, Nx, Ntime))
        for i in range(f2D_nbfiles):
            f2D_tmp = np.loadtxt(f2D_filelist[i])
            f2D_evol[:, :, i] = np.reshape(f2D_tmp[:, 2], (Nv, Nx))
        
        f2D_evol = np.array(f2D_evol)
        f2D_evol_true = True
    else:
        f2D_evol_true = False


    s['Phi1D_evol'] = Phi1D_evol
    s['Nx'] = Nx
    s['Ntime'] = Ntime
    s['xg'] = xg
    s['vg'] = vg
    s['dens_evol'] = dens_evol
    s['velo_evol'] = velo_evol
    s['temp_evol'] = temp_evol
    s['flux_evol'] = flux_evol
    s['f1Dx_evol'] = f1Dx_evol
    s['f1Dv_evol'] = f1Dv_evol
    s['Nv'] = Nv
    s['f2D_evol'] = f2D_evol
    s['f2D_evol_true'] = f2D_evol_true

    
    return s

if __name__ == "__main__":
    sl2d = read_VlasovPoiss()
