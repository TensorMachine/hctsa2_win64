% For compilation under Windows/Octave you need to install the package
% octavex.x-headers.
% Under Linux/Octave, compilation works out of the box.
%
% There is a big variety of platforms running Matlab/Octave around and the
% script is tested for:
%
% Linux 32bit      / Matlab R2008b, 7.7
% Linux 64bit      / Matlab R2010b, 7.11
% Windows 7  32bit / Matlab R14, 2004, 7.0
% Windows 7  64bit / Matlab R2011b, 7.13
% Windows XP 32bit / Matlab R2010b, 7.11
% OS X 10.7.5      / Matlab R2013b
% OS X 10.9        / Matlab R2013b
%
% Linux 32bit      / Octave 3.2.4
% OS X 10.9.1      / Octave 3.8.0 built with Homebrew
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2014-02-13.

OCTAVE = exist('OCTAVE_VERSION') ~= 0;        % check if we run Matlab or Octave

me = mfilename;                                            % what is my filename
mydir = which(me); mydir = mydir(1:end-2-numel(me));        % where am I located
cd(mydir)

fprintf('Compiling solve_chol.c ...\n')
if OCTAVE                                                               % Octave
  if ismac
                                           % Accelerate framework needed on OS X
    fprintf('mkoctfile --mex solve_chol.c "-Wl,-framework,Accelerate"\n')
    mkoctfile --mex solve_chol.c '-Wl,-framework,Accelerate'
  else
    fprintf('mkoctfile --mex solve_chol.c\n')
    mkoctfile --mex solve_chol.c
  end
  delete solve_chol.o
else                                                                    % Matlab
  if ispc                                                              % Windows
    % Link with -lmwlapack and let mex locate the import library.
    % (Carried across the GPML 3.5 -> 4.2 upgrade: 4.2 reintroduces the same
    % cc.Manufacturer path construction. Both versions build only solve_chol.c,
    % so this patched make.m applies unchanged.)
    %
    % The original built an explicit path from mex.getCompilerConfigurations:
    %   ospath = ['extern/lib/win64/', lower(cc.Manufacturer)]
    % For MinGW that yields .../win64/gnu/libmwlapack.lib, but MathWorks ships
    % the MinGW import libraries under .../win64/mingw64/. gcc was handed a
    % path that does not exist. The path also contains a space in
    % "Program Files", which is a second way for that construction to fail.
    try
      fprintf('mex -O -lmwlapack solve_chol.c -output solve_chol\n')
      eval('mex -O -lmwlapack solve_chol.c -output solve_chol')
    catch err
      fprintf(2,'  -lmwlapack failed (%s); trying an explicit library path.\n',err.message);
      try
        cc = mex.getCompilerConfigurations; cc = lower(cc(1,1).Manufacturer);
      catch
        cc = 'microsoft';
      end
      if strncmp(cc,'gnu',3), cc = 'mingw64'; end     % MathWorks' folder name
      if numel(strfind(computer,'64')), w = 'win64'; else, w = 'win32'; end
      lib = fullfile(matlabroot,'extern','lib',w,cc,'libmwlapack.lib');
      fprintf('mex -O solve_chol.c -output solve_chol "%s"\n',lib);
      mex('-O','solve_chol.c','-output','solve_chol',lib);   % function form quotes for us
    end
  else
    fprintf('mex -O -lmwlapack solve_chol.c\n')
    eval('mex -O -lmwlapack solve_chol.c')
  end
end
fprintf('Done!\n')