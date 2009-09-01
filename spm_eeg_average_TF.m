function D = spm_eeg_average_TF(S)
% Average each channel over trials or trial types, for time-frequency data
% FORMAT D = spm_eeg_average_TF(S)
%
% S         - optional input struct
% (optional) fields of S:
% S.D           - MEEG object or filename of M/EEG mat-file with epoched TF data
% S.circularise - flag that indicates whether average is straight (0) or
%                 vector (1) of phase angles.
% S.robust      - (optional) - use robust averaging (only for power)
%                 .savew  - save the weights in an additional dataset
%                 .bycondition - compute the weights by condition (1,
%                                default) or from all trials (0)
%                 .ks     - offset of the weighting function (default: 3)
%
% Output:
% D         - MEEG object (also written to disk)
%__________________________________________________________________________
%
% This function averages single trial time-frequency data within trial type.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_average_TF.m 3341 2009-09-01 14:23:49Z vladimir $

SVNrev = '$Rev: 3341 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG TF averaging'); spm('Pointer','Watch');

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

%-Check data type
%--------------------------------------------------------------------------
if ~strcmp(D.type, 'single')
    error('This function can only be applied to single trial data');
elseif ~strncmp(D.transformtype, 'TF', 2) % TF and TFphase
    error('This function can only be applied to time-frequency data.');
end

%-Get other parameters
%--------------------------------------------------------------------------
try
    circularise = S.circularise;
catch
    circularise = 0;
    if strcmp(D.transformtype, 'TFphase')
        circularise = spm_input('Average of phase angles?','+1', ...
            'straight|vector(PLV)',[0 1], 1);
    end
    S.circularise = circularise;
end


%-Configure robust averaging
%--------------------------------------------------------------------------

if strcmp(D.transformtype, 'TFphase')
    S.robust = 0;
end

if ~isfield(S, 'robust')
    robust = spm_input('Use robust averaging?','+1','yes|no',[1 0]);
    if robust
        S.robust = [];
    else
        S.robust = false;
    end
else
    robust = isa(S.robust, 'struct');
end


if robust
    if ~isfield(S.robust, 'savew')
        S.robust.savew =  spm_input('Save weights?','+1','yes|no',[1 0]);
    end
    savew = S.robust.savew;
    
    if ~isfield(S.robust, 'bycondition')
        S.robust.bycondition = spm_input('Compute weights by condition?','+1','yes|no',[1 0]);
    end
    
    bycondition = S.robust.bycondition;
    
    if ~isfield(S.robust, 'ks')
        S.robust.ks =  spm_input('Offset of the weighting function', '+1', 'r', '3', 1);
    end
    
    ks          = S.robust.ks;
    
end

%-Generate new MEEG object with new files
%--------------------------------------------------------------------------
Dnew = clone(D, ['m' fnamedat(D)], [D.nchannels D.nfrequencies D.nsamples D.nconditions]);

if robust && savew
    Dw = clone(D, ['W' fnamedat(D)], size(D));
end

cl   = D.condlist;

spm_progress_bar('Init', D.nchannels, 'Channels completed');

ni = zeros(1,D.nconditions);
for i = 1:D.nconditions
    w = pickconditions(D, deblank(cl{i}), 1)';
    ni(i) = length(w);
    if ni(i) == 0
        warning('%s: No trials for trial type %d', D.fname, cl{i});
    end
end

goodtrials  =  pickconditions(D, cl, 1);

for j = 1:D.nchannels
    if robust && ~bycondition
        [Y, W1] = spm_robust_average(D(j, :, :, goodtrials), 4, ks);
        if savew
            Dw(j, :, :, goodtrials) = W1;
        end
        W = zeros([1 D.nfrequencies D.nsamples D.ntrials]);
        W(1, :, :, goodtrials) = W1;
    end
    for i = 1:D.nconditions
        
        w = pickconditions(D, deblank(cl{i}), 1)';
        
        if isempty(w)
            continue;
        end
        
        %-Straight average
        %------------------------------------------------------------------
        if ~circularise
            Dnew(j, :, :, i) = mean(D(j, :, :, w), 4);
            
            %-Vector average (eg PLV for phase)
            %------------------------------------------------------------------
        else
            tmp = D(j, :, :, w);
            tmp = cos(tmp) + sqrt(-1)*sin(tmp);
            
            Dnew(j, :, :, i) = abs(mean(tmp,4)) ./ mean(abs(tmp),4);
        end
    end
    spm_progress_bar('Set', j);
end

spm_progress_bar('Clear');

Dnew = type(Dnew, 'evoked');

%-Update some header information
%--------------------------------------------------------------------------
Dnew = conditions(Dnew, [], cl);
Dnew = repl(Dnew, [], ni);

%-Display averaging statistics
%--------------------------------------------------------------------------
disp(sprintf('%s: Number of replications per contrast:', Dnew.fname));  %-#
s = [];
for i = 1:D.nconditions
    s = [s sprintf('average %s: %d trials', cl{i}, ni(i))];
    if i < D.nconditions
        s = [s sprintf(', ')];
    else
        s = [s '\n'];
    end
end
disp(sprintf(s));                                                       %-#

%-Save new evoked M/EEG dataset
%--------------------------------------------------------------------------
Dnew = Dnew.history(mfilename, S);
save(Dnew);

if robust && savew
    Dw = Dw.history(mfilename, S);
    save(Dw);
end

D = Dnew;

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG TF averaging: done'); spm('Pointer', 'Arrow');
