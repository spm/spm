function [x] = spm_load(f)
% function to load ascii file data as matrix
% FORMAT [x] = spm_load(f)
% f  - file {ascii file containing a regular array of numbers
% x  - corresponding data matrix
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_load.m 1143 2008-02-07 19:33:33Z spm $



%-Get a filename if none was passed
%-----------------------------------------------------------------------
x  = [];
if nargin == 0
    [f,p] = uigetfile({'*.mat';'*.txt';'*.dat'});
    try
        f     = fullfile(p,f);
    end
end

%-Load the data file into double precision matrix x
%-----------------------------------------------------------------------
try
    x = load(f,'-ascii');
    return
end
try
    x = load(f,'-mat');
    x = getdata(x);
end

if ~isnumeric(x), x = []; end

function x = getdata(s)
% get numberic data x from the fields of structure s
%--------------------------------------------------------------------------
x = [];
f = fieldnames(s);
for i = 1:length(f)
    x = s.(f{i});
    if isnumeric(x),return;         end
    if isstruct(x), x = getdata(x); end
end
