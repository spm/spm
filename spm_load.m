function [x] = spm_load(f)
% function to load ascii file data as matrix
% FORMAT [x] = spm_load(f)
% f  - file {ascii file containing a regular array of numbers
% x  - corresponding data matrix
%_______________________________________________________________________
% %E% Karl Friston, Andrew Holmes %W%


%-Get a filename if none was passed
%-----------------------------------------------------------------------
if nargin == 0
	[f,p] = uigetfile({'*.mat';'*.txt';'*.dat'});
    f     = fullfile(p,f);
end

%-Load the data file into double precision matrix x
%-----------------------------------------------------------------------
try
    x = load(f,'-ascii');
    return
end
try
    x = load(f,'-mat');
    while ~isnumeric(x)
        s = fieldnames(x);
        x = getfield(x,s{1});
    end
end
