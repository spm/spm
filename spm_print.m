function spm_print(F)
% SPM print function: Appends footnote & executes PRINTSTR
% FORMAT spm_print(F)
% F	- [Optional] Figure to print ('Tag' or figure number)
%	- defaults to the figure 'Tag'ged as 'Graphics'
%	- If no such figure found, uses CurrentFigure, if avaliable
%_______________________________________________________________________
%
% spm_print creates a footnote with details of the current
% session and evaluates the global string variable PRINTSTR
%
% This is just a gateway to spm_figure('Print',F)
%_______________________________________________________________________
% @(#)spm_print.m	1.3 96/04/25

if nargin==0
	spm_figure('Print')
else
	spm_figure('Print',F)
end
