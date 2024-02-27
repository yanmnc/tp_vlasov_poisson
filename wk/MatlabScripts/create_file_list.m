%------------------------------------------------------------
%  file : create_file_list.m
%  date : 2012-02-15
%   used for creating the list of the files 
%    <acronym>_<num>.extension
%------------------------------------------------------------
function [file_list,nb_files] = create_file_list(acronym,extension,syntax)
ifile     = 0;
filename  = create_file_name(acronym,0,extension,syntax);
nb_files  = 0;
file_list = [];
while (exist(filename)==2)
  ifile            = ifile + 1;
  nb_files         = nb_files+ 1;
  file_list{ifile} = filename;
  filename         = create_file_name(acronym,ifile,extension,syntax);
end
