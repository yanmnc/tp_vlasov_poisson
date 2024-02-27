%-----------------------------------------------------
%  file : read_VlasovPoiss.m
%  date : 2013/02/10
%
%   Save all the variables concerning the 2D Vlasov
%    Poisson simulation in a structure
%-----------------------------------------------------
function [sl2d] = read_VlasovPoiss()
     
%*** choice of the directory (if necessary) ***
Directory = input('Directory for data reading  [default = current] ? :','s');
if (~isempty(Directory))
  Prefix = [Directory,'/'];
else
  Prefix = '';
end

nb_s = 0;      %-> counter for the number of variables saved in the structure s

%*** Read the file 'param_simu.dat' ***
disp(['--> Reading of the param_simu.dat file']);
f_param      = load([Prefix,'param_simu.dat']);
name{nb_s+1} = 'kx0';
name{nb_s+2} = 'perturb';
name{nb_s+3} = 'v0';
name{nb_s+4} = 'T0';
name{nb_s+5} = 'epsilon';
name{nb_s+6} = 'vdiag_forf1Dx';
name{nb_s+7} = 'nu_coll';
name{nb_s+8} = 'Phi_ext0';
name{nb_s+9} = 'kx_ext';
name{nb_s+10} = 'omega_ext';
%nbparam = 10;
name{nb_s+11} = 'coef_OhmsLaw';
name{nb_s+12} = 'E0';
nbparam = 12;
for i=1:nbparam
  s.(name{i}) = f_param(i);
end  
nb_s = nb_s + nbparam;

%*** Read the file 'VlasovPoiss_res.dat' ***
disp(['--> Reading of the VlasovPoiss_res.dat file']);
f_pb2d       = load([Prefix,'VlasovPoiss_res.dat']);
name{nb_s+1} = 'time';
name{nb_s+2} = 'nbions';
name{nb_s+3} = 'Enkin';
name{nb_s+4} = 'Enpot';
name{nb_s+5} = 'entropy';
name{nb_s+6} = 'L1_norm';
name{nb_s+7} = 'L2_norm';   
nbpbeam     = 7;
for i=1:nbpbeam
  s.(name{nb_s+i}) = f_pb2d(:,i);
end  
nb_s = nb_s + nbpbeam;

%*** Read all the files 'Phi1D_<num>.dat' ***
disp(['--> Reading of the Phi1D_<num>.dat files']);
name{nb_s+1} = 'Nx';
name{nb_s+2} = 'Ntime';
name{nb_s+3} = 'xg';
name{nb_s+4} = 'Phi1D_evol';
Phi1D_acronym = [Prefix,'Phi1D'];
[Phi1D_filelist,Phi1D_nbfiles] = create_file_list(Phi1D_acronym,'.dat',0);
for i = 1:Phi1D_nbfiles,
  Phi1D_tmp       = load(Phi1D_filelist{i});
  Phi1D_evol(:,i) = Phi1D_tmp(:,2); 
end
s.('Phi1D_evol') = Phi1D_evol;
s.('Nx')         = size(Phi1D_evol,1);
s.('Ntime')      = size(Phi1D_evol,2);
s.('xg')         = Phi1D_tmp(:,1);
nb_s  = nb_s + 4;

%*** Read all the files 'dens_<num>.dat' ***
disp(['--> Reading of the dens_<num>.dat files']);
name{nb_s+1} = 'Nx';
name{nb_s+2} = 'Ntime';
name{nb_s+3} = 'xg';
name{nb_s+4} = 'dens_evol';
dens_acronym = [Prefix,'dens'];
[dens_filelist,dens_nbfiles] = create_file_list(dens_acronym,'.dat',0);
for i = 1:dens_nbfiles,
  dens_tmp       = load(dens_filelist{i});
  dens_evol(:,i) = dens_tmp(:,2); 
end
s.('dens_evol') = dens_evol;
s.('Nx')         = size(dens_evol,1);
s.('Ntime')      = size(dens_evol,2);
s.('xg')         = dens_tmp(:,1);
nb_s  = nb_s + 4;

%*** Read all the files 'velo_<num>.dat' ***
disp(['--> Reading of the velo_<num>.dat files']);
name{nb_s+1} = 'Nx';
name{nb_s+2} = 'Ntime';
name{nb_s+3} = 'xg';
name{nb_s+4} = 'velo_evol';
velo_acronym = [Prefix,'velo'];
[velo_filelist,velo_nbfiles] = create_file_list(velo_acronym,'.dat',0);
for i = 1:velo_nbfiles,
  velo_tmp       = load(velo_filelist{i});
  velo_evol(:,i) = velo_tmp(:,2); 
end
s.('velo_evol') = velo_evol;
s.('Nx')         = size(velo_evol,1);
s.('Ntime')      = size(velo_evol,2);
s.('xg')         = velo_tmp(:,1);
nb_s  = nb_s + 4;

%*** Read all the files 'temp_<num>.dat' ***
disp(['--> Reading of the temp_<num>.dat files']);
name{nb_s+1} = 'Nx';
name{nb_s+2} = 'Ntime';
name{nb_s+3} = 'xg';
name{nb_s+4} = 'temp_evol';
temp_acronym = [Prefix,'temp'];
[temp_filelist,temp_nbfiles] = create_file_list(temp_acronym,'.dat',0);
for i = 1:temp_nbfiles,
  temp_tmp       = load(temp_filelist{i});
  temp_evol(:,i) = temp_tmp(:,2); 
end
s.('temp_evol') = temp_evol;
s.('Nx')         = size(temp_evol,1);
s.('Ntime')      = size(temp_evol,2);
s.('xg')         = temp_tmp(:,1);
nb_s  = nb_s + 4;

%*** Read all the files 'flux_<num>.dat' ***
disp(['--> Reading of the flux_<num>.dat files']);
name{nb_s+1} = 'Nx';
name{nb_s+2} = 'Ntime';
name{nb_s+3} = 'xg';
name{nb_s+4} = 'flux_evol';
flux_acronym = [Prefix,'flux'];
[flux_filelist,flux_nbfiles] = create_file_list(flux_acronym,'.dat',0);
for i = 1:flux_nbfiles,
  flux_tmp       = load(flux_filelist{i});
  flux_evol(:,i) = flux_tmp(:,2); 
end
s.('flux_evol') = flux_evol;
s.('Nx')         = size(flux_evol,1);
s.('Ntime')      = size(flux_evol,2);
s.('xg')         = flux_tmp(:,1);
nb_s  = nb_s + 4;

%*** Read all the files 'f1Dx_<num>.dat' ***
disp(['--> Reading of the f1Dx_<num>.dat files']);
name{nb_s+1} = 'f1Dx_evol';
f1Dx_acronym = [Prefix,'f1Dx'];
[f1Dx_filelist,f1Dx_nbfiles] = create_file_list(f1Dx_acronym,'.dat',0);
for i = 1:f1Dx_nbfiles,
  f1Dx_tmp       = load(f1Dx_filelist{i});
  f1Dx_evol(:,i) = f1Dx_tmp(:,2); 
end
s.('f1Dx_evol') = f1Dx_evol;
nb_s = nb_s + 1;

%*** Read all the files 'f1Dv_<num>.dat' ***
disp(['--> Reading of the f1Dv_<num>.dat files']);
name{nb_s+1} = 'Nv';
name{nb_s+2} = 'vg';
name{nb_s+3} = 'f1Dv_evol';
f1Dv_acronym = [Prefix,'f1Dv'];
[f1Dv_filelist,f1Dv_nbfiles] = create_file_list(f1Dv_acronym,'.dat',0);
for i = 1:f1Dv_nbfiles,
  f1Dv_tmp       = load(f1Dv_filelist{i});
  f1Dv_evol(:,i) = f1Dv_tmp(:,2); 
end
s.('f1Dv_evol') = f1Dv_evol;
s.('Nv')        = size(f1Dv_evol,1);
s.('vg')        = f1Dv_tmp(:,1);
nb_s  = nb_s + 3;

%*** Read all the 2D files 'f2D_<num>.dat' ***
f2D_acronym = [Prefix,'f2D'];
[f2D_filelist,f2D_nbfiles] = create_file_list(f2D_acronym,'.dat',0);
name{nb_s+1} = 'f2D_evol';
if (f2D_nbfiles>1)
  Nx = s.('Nx');
  Nv = s.('Nv');
  disp(['--> Reading of the f2D files'])
  f2D_evol = [];
  for i = 1:f2D_nbfiles,
    f2D_tmp         = load(f2D_filelist{i});
    f2D_evol(:,:,i) = reshape(f2D_tmp(:,3),Nx,Nv); 
  end
  s.('f2D_evol') = f2D_evol;
  nb_s  = nb_s + 1;
  name{nb_s+1} = 'f2D_evol_true';
  s.('f2D_evol_true') = 1;
  nb_s  = nb_s + 1;
else
  name{nb_s+1} = 'f2D_evol_true';
  s.('f2D_evol_true') = 0;
  nb_s  = nb_s + 1;
end

%*** return the structure ***
if (nargout==0)
   for i=1:nb_s
     assignin('base',name{i},s.(name{i}));  
   end
else
   sl2d = s;
end
