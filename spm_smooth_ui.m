function spm_smooth_ui
% Smoothing or convolving
%___________________________________________________________________________
%
% Convolves image files with an isotropic (in real space) Gaussian kernel 
% of a specified width.
%
% Uses:
%
% As a preprocessing step to suppress noise and effects due to residual 
% differences in functional and gyral anatomy during inter-subject 
% averaging.
%
% Inputs
%
% *.img conforming to SPM data format (see 'Data')
%
% Outputs
%
% The smoothed images are written to the same subdirectories as the 
% original *.img and are prefixed with a 's' (i.e. s*.img)
%
%__________________________________________________________________________
% %W%	%E%

% Programmers Guide
% Batch system implemented on this routine. See spm_bch.man
% If inputs are modified in this routine, try to modify spm_bch.man
% and spm_bch_bchmat (if necessary) accordingly. 
% Calls to spm_input in this routine use the BCH gobal variable.  
%    BCH.bch_mat 
%    BCH.index0  = {'smooth',index_of_Analysis};
%_______________________________________________________________________

global BCH;

% get filenames and kernel width
%----------------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,'2.4');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Smooth');
spm_help('!ContextHelp','spm_smooth_ui.m');

s     = spm_input('smoothing {FWHM in mm}',1,...
                  'batch',{},'FWHMmm');
if isempty(BCH)
   P = spm_get(Inf,'.img','select scans');
else
   P = spm_input('batch',{},'files');
end
n     = size(P,1);

% implement the convolution
%---------------------------------------------------------------------------
spm('Pointer','Watch');
spm('FigName','Smooth: working',Finter,CmdLine);
spm_progress_bar('Init',n,'Smoothing','Volumes Complete');
for i = 1:n
	Q = deblank(P(i,:));
	[pth,nm,xt,vr] = fileparts(deblank(Q));
	U = fullfile(pth,['s' nm xt vr]);
	spm_smooth(Q,U,s);
	spm_progress_bar('Set',i);
end
spm_progress_bar('Clear',i);
spm('FigName','Smooth: done',Finter,CmdLine);
spm('Pointer');
