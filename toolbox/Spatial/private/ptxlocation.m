function ptx = ptxlocation(nam)
% Location of a PTX file used in GPU computations
% FORMAT ptx = ptxlocation(nam)
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


d = fileparts(mfilename('fullpath'));
ptx = fullfile(d,'..','lib',parallel.gpu.ptxext);

if nargin
    ptx = spm_file(nam,'path',ptx,'ext','.ptx');
end
