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
%     S.tf.Db          - MEEG object or filename of M/EEG mat-file to use
%                        for the baseline (if different from the input dataset).
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
% $Id: spm_eeg_tf_rescale.m 4316 2011-04-26 16:52:28Z vladimir $

SVNrev = '$Rev: 4316 $';

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
    str  = {'LogR','Diff', 'Rel', 'Log', 'Sqrt', 'Zscore'};
    S.tf.method = spm_input('Rescale method','+1','m',str,char(str),1);
end

Din  = spm_eeg_load(D);
tims = time(Din);

Nf   = length(frequencies(Din));
D    = clone(Din, ['r' Din.fnamedat], [Din.nchannels Nf Din.nsamples Din.ntrials]);

switch lower(S.tf.method)
    
    case {'logr','diff', 'rel', 'zscore'}
        try
            S.tf.Sbaseline;
        catch            
            if spm_input('Baseline dataset','+1','b',{'Same|Different'},[0 1],0)
                [Db, sts] = spm_select(1, 'mat', 'Select baseline M/EEG mat file');
                if ~sts, return; end
                S.tf.Db = Db;
            else
                S.tf.Db = [];
            end
            
            tmp_base = spm_input('Start and stop of baseline [ms]', '+1', 'i', '', 2);
            S.tf.Sbaseline = tmp_base/1000;
        end
        
        if isfield(S.tf, 'Db') && ~isempty(S.tf.Db)
            Db = spm_eeg_load(S.tf.Db);
        else
            Db = Din;
        end
        
        if any(abs(Din.frequencies-Db.frequencies)>0.1) || ~isequal(Db.chanlabels, Din.chanlabels) ||...
                (Db.ntrials>1 && (Db.ntrials~=Din.ntrials))
            error('The input dataset and the baseline dataset should have the same frequencies, channels and trial numbers');
        end        
       
        for c=1:D.ntrials
            inds=find(tims>=S.tf.Sbaseline(1) & tims<=S.tf.Sbaseline(2));
            x=spm_squeeze(Din(:,:,:,c), 4);
            if Db.ntrials > 1
                xbase=spm_squeeze(Db(:,:,:,c), 4);
            else
                xbase=spm_squeeze(Db(:,:,:,1), 4);
            end            
            switch lower(S.tf.method)
                case 'logr'
                    xbase=mean(log10(xbase(:,:,inds)),3);
                    D(:,:,:,c)= 10*(log10(x) - repmat(xbase,[1 1 D.nsamples 1]));
                    D = units(D, [], 'dB');
                case 'diff'
                    xbase=mean(xbase(:,:,inds),3);
                    D(:,:,:,c)= (x - repmat(xbase,[1 1 D.nsamples 1]));
                case 'zscore'
                    stdev = std(xbase(:,:,inds), [], 3);
                    xbase= mean(xbase(:,:,inds),3);                    
                    D(:,:,:,c)= (x - repmat(xbase,[1 1 D.nsamples 1]))./repmat(stdev,[1 1 D.nsamples 1]);
                case 'rel'
                    xbase=mean(xbase(:,:,inds),3);
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
