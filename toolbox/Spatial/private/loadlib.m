function loadlib(nam)
% Load a shared library into MATLAB
% FORMAT loadlib(nam)
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


if ~nargin, nam = 'pushpull'; end

if libisloaded(nam), return; end

d = fileparts(mfilename('fullpath'));
mext  = mexext;
soext = spm_platform('soext');
libname = fullfile(d,'..','lib',['lib' mext(4:end)],[nam '.' soext]);
hdrname = fullfile(d,'..','src',  [nam '.h']);

[notfound,warnings] = loadlibrary(libname,hdrname);

if ~isempty(notfound)
    fprintf('\n --- Displaying debugging information ---\n');
    disp(notfound)
    disp(warnings)
    fprintf('\n ---                                  ---\n');
end
