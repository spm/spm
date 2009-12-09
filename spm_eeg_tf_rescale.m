function [D] = spm_eeg_tf_rescale(S)
% Rescale (avg) spectrogram with nonlinear and/or difference operator
% FORMAT [D] = spm_eeg_tf_rescale(S)
%
% S                    - input structure (optional)
% (optional) fields of S:
%   S.D                - MEEG object or filename of M/EEG mat-file
%   S.tf               - structure with (optional) fields:
%     S.tf.method      - 'LogR', 'Diff', 'Rel', 'Log', 'Sqrt'
%     S.tf.Sbaseline   - 2-element vector: start and stop of baseline 
%                        (need to specify this for LogR and Diff)
% 
% D                    - MEEG object with rescaled power data (also
%                        written to disk with prefix r)
%
% For 'Log' and 'Sqrt', these functions are applied to spectrogram 
% For 'LogR', 'Rel' and 'Diff' this function computes power in the baseline
% p_b and outputs (i) p-p_b for 'Diff' (ii) 100*(p-p_b)/p_b for 'Rel' 
%                 (iii) log (p/p_b) for 'LogR'
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_eeg_tf_rescale.m 3622 2009-12-09 09:36:39Z vladimir $

SVNrev = '$Rev: 3622 $';

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
    str  = {'LogR','Diff', 'Rel', 'Log','Sqrt'};
    S.tf.method = spm_input('Rescale method','+1','b',str,[],1);
end

Din  = spm_eeg_load(D);
tims = time(Din);

Nf   = length(frequencies(Din));
D    = clone(Din, ['r' Din.fnamedat], [Din.nchannels Nf Din.nsamples Din.ntrials]);

switch lower(S.tf.method)
    
    case {'logr','diff', 'rel'}
        try
            S.tf.Sbaseline;
        catch
            tmp_base = spm_input('Start and stop of baseline [ms]', '+1', 'i', '', 2);
            S.tf.Sbaseline = tmp_base/1000;
        end
        for c=1:D.ntrials
            inds=find(tims>=S.tf.Sbaseline(1) & tims<=S.tf.Sbaseline(2));
            % reshape instead of squeeze as there may be other singleton
            % dimensions
            x=reshape(Din(:,:,:,c), size(Din, [1:3]));
            xbase=mean(x(:,:,inds),3);
            switch lower(S.tf.method)
                case 'logr'
                    x=log10(x);
                    xbase=mean(x(:,:,inds),3);
                    D(:,:,:,c)= 10*(x - repmat(xbase,[1 1 D.nsamples 1]));
                    D = units(D, [], 'dB');
                case 'diff'
                    D(:,:,:,c)= (x - repmat(xbase,[1 1 D.nsamples 1]));
                case 'rel'
                    D(:,:,:,c)= 100*((x./repmat(xbase,[1 1 D.nsamples 1]) - 1));
                    D = units(D, [], '%');
            end
        end
        
    case 'log'
        for c=1:D.ntrials
            D(:,:,:,c) = log(Din(:,:,:,c));
        end
        
    case 'sqrt'
        for c=1:D.ntrials
            D(:,:,:,c) = sqrt(Din(:,:,:,c));
        end
        
    otherwise
        error('Unknown rescaling method.');
end

% Save
D = D.history(mfilename, S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG Time Frequency Rescale: done'); 
spm('Pointer','Arrow');
