function [D] = spm_eeg_tf_rescale(S)
% Rescale (avg) spectrogram with nonlinear and/or difference operator
% FORMAT [D] = spm_eeg_tf_rescale(S)
%
% S                    - input structure (optional)
% (optional) fields of S:
%   S.D                - MEEG object or filename of M/EEG mat-file
%   S.tf               - structure with (optional) fields:
%     S.tf.method      - 'LogR', 'Diff', 'Log', 'Sqrt'
%     S.tf.Sbaseline   - 2-element vector: start and stop of baseline 
%                        (need to specify this for LogR and Diff)
% 
% D                    - MEEG object with rescaled power data (also
%                        written to disk with prefix r)
%
% For 'Log' and 'Sqrt', these functions are applied to spectrogram 
% For 'LogR' and 'Diff' this function computes power in the baseline
% p_b and outputs (i) p-p_b for 'Diff' or (ii) log (p/p_b) for 'LogR'
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_eeg_tf_rescale.m 3089 2009-04-29 17:24:39Z will $

SVNrev = '$Rev: 3089 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG Time-Frequency Rescale'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

try
    S.tf.method;
catch
    str  = {'LogR','Diff','Log','Sqrt'};
    S.tf.method = spm_input('Rescale method','+1','b',str,[],1);
end

Din = spm_eeg_load(D);
tims=time(Din);

Nf=length(frequencies(Din));
D = clone(Din, ['r' Din.fnamedat], [Din.nchannels Nf Din.nsamples Din.ntrials]);

switch lower(S.tf.method),
    case {'logr','diff'}
        try
            S.tf.Sbaseline;
        catch
            tmp_base = spm_input('Start and stop of baseline [ms]', '+1', 'i', '', 2);
            S.tf.Sbaseline = tmp_base/1000;
        end
        for c=1:D.ntrials,
            inds=find(tims>S.tf.Sbaseline(1) & tims<S.tf.Sbaseline(2));
            x=squeeze(Din(:,:,:,c));
            xbase=mean(x(:,:,inds),3);
            if strcmp(lower(S.tf.method),'logr')
                x=log(x);
                xbase=log(xbase);
            end
            D(:,:,:,c)= x - repmat(xbase,[1 1 D.nsamples 1]);
        end
    case 'log',
        for c=1:D.ntrials,
            fx=log(squeeze(Din(:,:,:,c)));
            D(:,:,:,c)=fx;
        end
    case 'sqrt',
        for c=1:D.ntrials,
            fx=sqrt(squeeze(Din(:,:,:,c)));
            D(:,:,:,c)=fx;
        end
    otherwise
        disp('Error in spm_eeg_tf_rescale: unknown method');
        return;
end

% Save
D = D.history(mfilename, S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG Time Frequency Rescale: done'); 
spm('Pointer','Arrow');
