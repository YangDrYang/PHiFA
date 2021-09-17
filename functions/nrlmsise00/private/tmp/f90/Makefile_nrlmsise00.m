%ALL: mex

clc;
path = fileparts(mfilename('fullpath'));
cd(path);
%mex:
mex_file = 'nrlmsise00';
clear *.mex* 
delete *.mex*

mex('-v','-compatibleArrayDims','-O',...
  'nrlmsise00_constants.f90', 'm_nrlmsise00.f90',...
  [mex_file,'_mex.F'],...
  '-output',[mex_file,'_mex']);
% mex -v -compatibleArrayDims -g ...
%   FFLAGS='$FFLAGS -ffpe-trap=invalid,zero,overflow ' -g ...
%   f77/m_nrlmsise00.f  f77/nrlmsise00_mex.F...
%   -output  nrlmsise00_mex
  
%clean:
delete *.mod  *.*~  *.o

%install:
mexfile = dir('*.mex*');
if ~isempty(mexfile)
  disp([mex_file,' compile successfully!']);
  % move to the upper dir
  mexfilename = mexfile(1).name;
  delete(fullfile('..',mexfilename));
  movefile(mexfilename, '../');
  delete('../../nrlmsise00.m');
  copyfile('nrlmsise00.m','./../../','f');
else
  disp([mex_file,' compile failed!']);
end