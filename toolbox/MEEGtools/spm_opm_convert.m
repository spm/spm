function Dout = spm_opm_convert(array,fnamedat,fs,scale)
% Convert array into SPM MEEG object
% FORMAT Dout = spm_opm_convert(array,fnamedat,fs,scale)
%
% array       - numeric Array of 2,3,4 dimensions(channels,time,trials)
% fnamedat    - string specifing output path of object(include extension .dat)
% fs          - sampling frequency
% scale       - scale factor to convert to fT [default: 1]
%__________________________________________________________________________
% Copyright (C) 2017-2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_opm_convert.m 7414 2018-09-07 11:00:29Z spm $

% determine output filename
[a, b] = fileparts(fnamedat);
outMat = fullfile(a,[b,'.mat']);

% if scale is not supplied set a default  of 1 
if nargin < 4
    scale = 1;
end
    
array = array.*scale;

% find number of dimensions and decide what to do based on result
dim = size(array);
L   = length(dim);

if L > 4
    % throw error for having unsupported number of dimensions
    error('Array should not have more than 4 dimensions.');
elseif L == 4
    % create MEG object with appropriate size
    Dout = meeg(dim(1),dim(2),dim(3),dim(4));
    
    % make it a blank object
    Dout= blank(Dout,fnamedat);
    
    % set the smaple rate 
    Dout = Dout.fsample(fs);
    
    % fill in data with supplied array
    Dout(1:dim(1),1:dim(2),1:dim(3),1:dim(4)) = array;
elseif L == 3
    % same with less dimensions
    Dout = meeg(dim(1),dim(2),dim(3));
    Dout= blank(Dout,fnamedat);
    Dout = Dout.fsample(fs);
    Dout(1:dim(1),1:dim(2),1:dim(3)) = array;
    
elseif L == 2 
    %same with less dimensions
    Dout = meeg(dim(1),dim(2),1);
    Dout= blank(Dout,fnamedat);
    Dout = Dout.fsample(fs);
    Dout(1:dim(1),1:dim(2),1) = array;
else
    % throw exception if I really don't know what to do
    error('Array must have between 2 and 4 dimensions.');
end

% set the filename and save
Dout = fname(Dout,outMat);
Dout.save;
