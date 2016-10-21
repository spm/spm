function D = spm_eeg_erp_correction(S)
% Applies corrections to ERPs or single trials as in DCM-ERP
% This can be used to make a sensor level analysis or source reconstruction
% consistent with DCM.
% FORMAT D = spm_eeg_erp_correction(S)
%
% S        - optional input struct
% (optional) fields of S:
% S.D        - MEEG object or filename of M/EEG mat-file with epoched data
% S.detrend  - detrending order (0 for no detrending)
% S.hanning  - apply Hanning window (true or false)
%
% Output:
% D        - MEEG object (also written on disk)
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Melanie Boly
% $Id: spm_eeg_erp_correction.m 6907 2016-10-21 09:41:59Z vladimir $

SVNrev = '$Rev: 6907 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Correct ERPs'); spm('Pointer','Arrow');

if nargin == 0
    S = [];
end

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);


if isequal(D.type, 'continuous')
    error('Corrections can only be applied to epoched or averaged datasets.');
end


%-Get parameters
%--------------------------------------------------------------------------
if ~isfield(S, 'detrend')
    S.detrend = spm_input('Detrend order', '+1', 'n', '2', 1);
end

if ~isfield(S, 'hanning')
    S.hanning = spm_input('Apply Hanning?','+1','yes|no',[1 0], 0);
end

Ns = D.nsamples;

%-Confounds - DCT:
%--------------------------------------------------------------------------
if S.detrend == 0
    X0 = sparse(Ns, 1);
else
    X0 = spm_dctmtx(Ns, S.detrend);
end
R      = speye(Ns) - X0*X0';

%-Hanning
%--------------------------------------------------------------------------
if S.hanning
    R  = R*diag(spm_hanning(Ns))*R;
end

Dnew = clone(D, ['C' fname(D)]);


spm_progress_bar('Init', D.ntrials, 'Trials filtered'); drawnow;
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
else Ibar = [1:D.ntrials]; end

spm('Pointer','Watch');

%-Adjust data
%--------------------------------------------------------------------------
for i = 1:D.ntrials
    Dnew(Dnew.indchantype('MEEG', 'GOOD'),:,i) = (R*spm_squeeze(D(D.indchantype('MEEG', 'GOOD'),:,i), 3)')';
    
    if ismember(i, Ibar), spm_progress_bar('Set', i); end
end

spm_progress_bar('Clear');

D = Dnew;

D = D.history(mfilename, S);

save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','Correct ERPs: done'); spm('Pointer','Arrow');
