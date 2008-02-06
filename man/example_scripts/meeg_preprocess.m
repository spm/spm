% preprocess 
% This is an example batch script for preprocessing multiple subjects' EEG data. 
% The data is not part of the SPM5 distribution. This script is intended as
% a template for similar batch scripts.

clear all

% Enter Directories and other parameters here
%**************************************
Pbase = 'D:\forpeople\marta\';
Pdata = fullfile(Pbase, 'alldata');
%**************************************

cd(Pdata);


Dpath = spm_select(Nsub, 'dir', 'Select subject directories');
Nsub = size(Dpath, 1);

for i = 1:Nsub
	
	clear S
    
    S.Fchannels = fullfile(spm('dir'), 'EEGtemplates', 'bdf_setup.mat');
    S.Fdata = deblank(Dpath(i, :));
	D = spm_eeg_rdata_bdf(S);


    % epoch
    clear S
    S.events.start = -100;
    S.events.stop = 400;
    S.events.types = [1 3];
    S.events.Inewlist = 0;
    S.D = D.fname;

    D = spm_eeg_epochs(S);

    % downsample
    clear S
    S.Radc_new = 200;
    S.D = D.fname;

    D = spm_eeg_downsample(S);


    % bandpass filter
    clear S
    S.filter.type = 'butterworth';
    S.filter.PHz = [0.5 40];
    S.filter.band = 'bandpass';
    S.D = D.fname;

    D = spm_eeg_filter(S);


    % artifact detection
    clear S
    S.artefact.weighted = 0;
    S.thresholds.Check_Threshold = 1;
    S.thresholds.threshold = 100;
    S.thresholds.External_list = 0;
    S.D = D.fname;

    D = spm_eeg_artefact(S);


    % average
    clear S
    S.D = D.fname;

    D = spm_eeg_average(S);


end

for i = 1:Nsub
    clear S
    S.Fname = fullfile(fileparts(fname{i}), ['mafde_' spm_str_manip(fname{i}, 'rt')]);
    S.n = 64;
    S.interpolate_bad = 1;
    spm_eeg_convertmat2ana(S);    
end
