function spm_fmri_realign_ui
% user interface routine for spm_fmri_realign
% FORMAT spm_fmri_realign_ui
%____________________________________________________________________________
%
% spm_fmri_realign_ui prompts for a list of filenames representing a series
% of fMRI scans that are to be realigned and invokes spm_fmri_realign.m
% The first image for each subject (or set of realignments) is taken to be
% the reference image (to which the remaining images are realigned)
%
%__________________________________________________________________________
% %W% Karl Friston %E%


%----------------------------------------------------------------------------
figure(2); clf
set(2,'Name','Realignment and adjustment of fMRI data')
FILTER = [];


% cycle through subjects and images
%----------------------------------------------------------------------------
n = spm_input('number of subjects',1);
for i = 1:n
	q = spm_input(['number of scans: subject ' num2str(i)],2);
	P = spm_get(q,'.img',['select scans for subject ' num2str(i)]);
	eval(['P' num2str(i) ' = P;']);
end


% implement spatial realignment
%----------------------------------------------------------------------------
set(2,'Name','executing','Pointer','Watch'); drawnow
for i = 1:n
	eval(['P = P' num2str(i) ';'])
	spm_fmri_realign(P);
end
figure(2); clf; set(2,'Name','','Pointer','Arrow'); 
