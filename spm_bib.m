function spm_bib
% creates a figure and displays spm_bib.man
% FORMAT spm_bib
%___________________________________________________________________________
%
% spm_bib simply displays release notes and a short SPM bibliography
% in a specially configured window.
%
% It is called by spm.m
%
%__________________________________________________________________________
% %W% %E%


%---------------------------------------------------------------------------
close all

% position figures and set PaperPosition for printing
%----------------------------------------------------------------------------
S    = get(0,'ScreenSize');
A    = diag([S(3)/1152 S(4)/900 S(3)/1152 S(4)/900]);
S3   = [276 008 600 865]*A;

whitebg(0,'w')

figure('Name','Release notes and Bibliogrpahy',...
	'NumberTitle','off','Position',S3,'Resize','off',...
	'PaperPosition',[.75 1.5 7 9.5],'Tag','Graphics')

set(gcf,'Pointer','Watch')
spm_help_disp('spm_bib.man')
set(gcf,'Pointer','Arrow')
