function [x] = spm_load(f)
% function to load ascii file data
% FORMAT [x] = spm_load(f)
% f  - file {ascii file containing a regular array of numbers
% x  - corresponding data matrix
%___________________________________________________________________________
% %E% Karl Friston, Andrew Holmes %W%

if nargin==0, f=spm_get(1,'*','Select ASCII data file'); end
eval(['load ' f ' -ascii; x = ' spm_str_manip(spm_str_manip(f,'t'),'r') ';']);
