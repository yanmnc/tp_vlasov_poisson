# Code for solving2D Vlasov-Poisson system

## Description of the code architecture

The base directory of the code repository contains the following sub-directories: 
 - `src`: contains the Fortran90 source files of the code;
 - `wk`: contains submission scripts to launch the code;
 - `post_process`: contains scripts to analyze the simulation results.

## Compilation
Use the following command to compile the `VlasovPoisson.exe` executable: 
```
cd src 
make clean
make
```

**Remark:** the compiler choice and the compiling options are done in the `Makefile` file. This file does not need to be changed as it already contains all required commands to compile the code properly. 


## Running the simulation
The simulation parameters are defined in submission scripts. Examples of such submission scripts are `wk/run_landau` and `wk/run_pbeam`. These scripts contain:
- A command for creating automatically a directory used for storing the results. For instance, if 
   the submission script is called `run_landau`, the corresponding result directory will be `DRUN_LANDAU`;
- The definition of the parameters of the simulation;
- The execution command.

Use the following command to run a submission script:
```
cd wk
./run_landau   # results will be saved in wk/DRUN_LANDAU/ 
```

or 
```
cd wk
./run_bot   # results will be saved in wk/DRUN_BOT/ 
```

**Remark.** The idea is to create a script per parameter set in order to store each simulation in a different directory


## Output results
The output files of the simulation are the following:  
1. `VlasovPoiss_res.out`: screen output. If you want to follow the time evolution of the simulation, use the command 
```
tail -f <RESULT_DIR>/VlasovPoiss_res.out
```
2. `VlasovPoiss_res.dat` : ASCII file containing the time evolution of the following quantities:
    - number of ions;
    - kinetic energy;
    - potential energy;
    - entropy;
    - L1 norm;  
    - L2 norm;
3. `Phi1D_<num_diag>.dat`: profile of the electrostatic potential, i.e Phi(x), 
                            at diagnostic `<num_diag>`;
4. `f1Dx_<num_diag>.dat`: x-plane of the distribution function, i.e f(x,v=fixed);
5. `f1Dv_<num_diag>.dat`: v-plane of the distribution function, i.e f(x=fixed,v);
6. `f2D_<num_diag>.dat`: (x,v) cross-section of the distribution function, i.e f(x,v).


## Post-process tools

### Matlab scripts
These scripts are located in the `post_process/MatlabScripts/` directory. The existing scrits are:
- `read_VlasovPoiss.m`: to read the Vlasov Poisson results;
- `PostProcess_Landau.m`: (or `PostProcess_BoumpOnTail.m`) : examples of plots;
- `Derivee1.m`: for 1D derivative computation;
- `Derivee2.m`: for 2D derivative computation;
- `Fourier1.m`: for Fast Fourier Transform in 1D;
- `Fourier2.m`: for Fast Fourier Transform in 2D.

To use one of the above scripts: 
```
module avail # (to see which modules are available for download)
module load matlab/XX # e.g. matlab/2013b
matlab #(to launch Matlab)
addpath('post_process/MatlabScripts/')
read_VlasovPoiss()           
PostProcess_Landau (or PostProcess_BoumpOnTail)
```

**Remarks:** 
- If you use the MATLAB command `read_VlasovPoiss()`, all the variables will be global variables. 
- If you use the MATLAB command `resu = read_VlasovPoiss()`, all the variables will be stored in the data structure `resu`. Hence, to access to the potential energy for instance, you have to write: `resu.Enpot` (instead of directly Enpot). This solution is useful if you want to analyse several simulations at the same time. 
