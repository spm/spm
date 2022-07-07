function [pth,nam,ext,num] = spm_fileparts(fname)
% Like fileparts, but separates off a comma separated list at the end
% FORMAT [pth,nam,ext,num] = spm_fileparts(fname)
% fname  - original filename
%
% pth    - path
% nam    - filename
% ext    - extension
% num    - comma separated list of values
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


num = '';
if ~ispc, fname = strrep(fname,'\',filesep); end
[pth,nam,ext] = fileparts(fname);
ind = find(ext==',');
if ~isempty(ind)
    num = ext(ind(1):end);
    ext = ext(1:(ind(1)-1));
end
