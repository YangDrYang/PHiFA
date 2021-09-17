%ALL: mex

clc;
path = fileparts(mfilename('fullpath'));
cd(path);
%mex:
mex_file = 'nrlmsise00';
clear *.mex* 
delete *.mex*

% mex('-v','-compatibleArrayDims','-g',...
%   [mex_file,'_mex.F'],...
%   'utils_constants.f90','utils_spline.f90', 'physics_constants.f90',...
%   'physics_msis.f90','nrlmsise00_constants.f90', 'm_nrlmsise00.f90',...
%   '-output',[mex_file,'_mex']);

mex('-v','-compatibleArrayDims','-O',...
  'm_nrlmsise00.f','nrlmsise00_mex.F',...
  '-output',[mex_file,'_mex']);

% mex('-v', '-compatibleArrayDims', '-g', ...
%   ['-I"',path,'"'],...
%   'CFLAGS="\$CFLAGS -Wcomment "',...%-std=c99
%   'nrlmsise-00.c', 'nrlmsise00_mex.c',...
%   '-output', 'nrlmsise00_mex');
  
%clean:
delete *.mod  *.*~  *.o

%install:
mexfile = dir('*.mex*');
if ~isempty(mexfile)
  disp([mex_file,' compile successfully!']);
  % move to the upper dir
  mexfilename = mexfile(1).name;
  movefile(mexfilename, './../');
  copyfile('nrlmsise00.m','./../../','f');
else
  disp([mex_file,' compile failed!']);
end