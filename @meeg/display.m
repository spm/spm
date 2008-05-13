function display(this)
% Method for displaying information about an meeg object
% FORMAT display(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: display.m 1621 2008-05-13 18:43:32Z vladimir $

disp('SPM M/EEG data object');
disp(['Type: ' type(this)]);
disp([num2str(nconditions(this)), ' conditions']);
disp([num2str(nchannels(this)), ' channels']);
disp([num2str(nsamples(this)), ' samples/trial']);
disp([num2str(ntrials(this)), ' trials']);
disp(['Sampling frequency: ' num2str(fsample(this)) ' Hz']);
disp(['Loaded from file ' fullfile(this.path, this.fname)]);
disp(' ');
disp('Use the syntax D(channels, samples, trials) to access the data');
disp('Type "methods(''meeg'')" for the list of methods performing other operations with the object');