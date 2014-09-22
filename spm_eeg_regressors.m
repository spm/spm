function regfile = spm_eeg_regressors(S)
% Prepare regressors for GLM analysis of M/EEG data
% FORMAT regfile = spm_eeg_regressors(S)
%
% S                     - input structure
%
% fields of S:
%   S.D                 - MEEG object or filename of M/EEG mat-file 
%                         for which the regressors should be prepared
%
% Output:
% regfile              - path to mat file in which the regressors are saved
%__________________________________________________________________________
% This is a modular function for which plugins can be developed implementing
% specific regressor creation cases
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_regressors.m 6186 2014-09-22 11:31:11Z vladimir $

SVNrev = '$Rev: 6186 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG Regressors'); spm('Pointer','Watch');

D = spm_eeg_load(S.D);

R = [];
names = {};

for i = 1:numel(S.methods)
    fun  = char(fieldnames(S.methods{i}));
    S1   = S.methods{i}.(fun);
    S1.D = D;
    S1.summarise = S.summarise;
    res =  feval(['spm_eeg_regressors_' fun], S1);
    
    R     =     [R res.R];
    names = [names res.names];
end


[outpath, outname] = fileparts(S.outfile);
if isempty(outpath)
    outpath = D.path;
end

regfile = fullfile(outpath, [outname '.mat']);

save(regfile, 'R', 'names', spm_get_defaults('mat.format'));

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG Regressors: done'); spm('Pointer','Arrow');
