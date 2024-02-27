%------------------------------------------------------------
%  file : create_file_name.m
%  date : 19/01/2006
%   used for creating the file name according 
%   to the number of the file
%  modification : 
%------------------------------------------------------------
function [file_name] = create_file_name(acronym,ifile,extension,syntax)
if (syntax==1)
  if (ifile<10)
    file_name = [acronym,'_d0000',num2str(ifile),extension];
  else 
    if (ifile<100) 
      file_name = [acronym,'_d000',num2str(ifile),extension];
    else
      if (ifile<1000)
	file_name = [acronym,'_d00',num2str(ifile),extension];
      else
	if (ifile<10000)
	file_name = [acronym,'_d0',num2str(ifile),extension];
	else
	file_name = [acronym,'_d',num2str(ifile),extension];
	end
      end
    end
  end
else
  if (ifile<10)
    file_name = [acronym,'_000',num2str(ifile),extension];
  else 
    if (ifile<100) 
      file_name = [acronym,'_00',num2str(ifile),extension];
    else
      if (ifile<1000)
	file_name = [acronym,'_0',num2str(ifile),extension];
      else
	file_name = [acronym,'_',num2str(ifile),extension];
      end
    end
  end
end
