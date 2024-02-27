!****************************************************************
! file : README.txt
! date : 2013-02-10
!
! README for Vlasov Poisson 2D code
!****************************************************************

Let us call DIR_VP the main directory for the Vlasov Poisson code.
(For instance DIR_VP=VlasovPoisson)
Then the main directory contains several sub-directories
 - src/            : for the Fortran90 source files of the code
 - wk/             : for the submission scripts
 - MatlabScripts/  : for MATLAB scripts --> examples to analyse the results


!--------------------------------------
! 1) COMPILING 
!--------------------------------------
Creation of the executable: 'VlasovPoisson.exe'
   > cd src 
   > make clean
   > make

Rk: 
1) The compiler choice and the compiling options are done in the 'Makefile' file.
2) Compilation has already been done. Unless you modify the source files, there is
no need to compile again


!--------------------------------------
! 2) PROGRAM RUNNING
!--------------------------------------
The simulation parameters are defined in a submission script (see for instance 
'run_landau' or 'run_pbeam'. These files contain:
  a) The automatic creation of a directory for the result storage. For instance if 
     the submission script is called 'run_landau', the corresponding result 
     directory will be DRUN_LANDAU.
  b) The definition of the parameters of the simulation
  c) The execution command
To use one of these scripts:
    > cd wk
    > ./run_landau &  ( results will be saved in wk/DRUN_LANDAU/ )
or  > ./run_bot  &  ( results will be saved in wk/DRUN_BOT/  )

Rk: The idea is to create a script per parameter set in order to store each 
    simulation in a different directory


!--------------------------------------
! 3) OUTPUT RESULTS
!--------------------------------------
The output files of the simulation are the following:
    a) VlasovPoiss_res.out : screen output. 
        If you want to follow the time evolution of the simulation, use the command:
	  tail -f <RESULT_DIR>/VlasovPoiss_res.out
    b) VlasovPoiss_res.dat : ASCII file containing the time evolution of the following 
       quantities:
        . ion number
        . kinetic energy
        . potential energy
        . entropy
        . L1 norm  
        . L2 norm
    c) Phi1D_<num_diag>.dat : profile of the electrostatic potential, i.e Phi(x), 
                                at diagnostic <num_diag>
    d) f1Dx_<num_diag>.dat  : x-plane of the distribution function, i.e f(x,v=fixed)
    e) f1Dv_<num_diag>.dat  : v-plane of the distribution function, i.e f(x=fixed,v)
    f) f2D_<num_diag>.dat   : (x,v) cross-section of the distribution function,
                                i.e f(x,v)


!--------------------------------------
! 4) TOOLS FOR RESULT ANALYSIS
!--------------------------------------
MATLAB scripts (defined in 'wk/MatlabScripts/' directory):
    - read_VlasovPoiss.m : to read the Vlasov Poisson results
    - PostProcess_Landau.m (or PostProcess_BoumpOnTail.m) : examples of plots 
    - Derivee1.m         : for 1D derivative computation
    - Derivee2.m         : for 2D derivative computation
    - Fourier1.m         : for Fast Fourier Transform in 1D
    - Fourier2.m         : for Fast Fourier Transform in 2D

   > module avail (to see which modules are available for download)
   > module load matlab/XX (e.g. matlab/2013b)
   > matlab                       (to launch Matlab)
   >> addpath('~/VlasovPoisson/wk/MatlabScripts/')
   >> read_VlasovPoiss()           
   >> PostProcess_Landau (or PostProcess_BoumpOnTail)

Remark: 
- If you use the MATLAB command:
    >> read_VlasovPoiss() 
  All the variables will be global variables. 

- If you use the MATLAB command: 
    >> resu = read_VlasovPoiss() 
  All the variables are stored in the structure 'resu'.
  For instance, to access to the potential energy, you have to write:
    >> resu.Enpot     (instead of directly Enpot)
  This solution is useful if you want to analyse several simulations 
  at the same time. 
