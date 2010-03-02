function Sout = spm_eeg_compatibility(S, caller)
% Keep compatibility for scripts by translating from old to new configuration
% FORMAT Sout = spm_eeg_compatibility(S, caller)
%
% Input:
%  S      - configuration struct
%  caller - calling function
%
% Output:
%  Sout   - updated configuration struct
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_compatibility.m 3742 2010-03-02 15:15:43Z vladimir $

switch caller
    case 'spm_eeg_tf'
        if isfield(S, 'tf')
            warning('The provided configuration of spm_eeg_tf is obsolete. Trying to update.');
            
            D = spm_eeg_load(S.D);
            
            Sout = [];
            
            Sout.D = S.D; 
            
            Sout.frequencies = S.tf.frequencies;                       
            Sout.channels    = D.chanlabels(S.tf.channels);
           
            if isfield(S.tf, 'pow')
                warning('''pow'' option is deprecated'); 
            end
            
            if isfield(S.tf, 'collchans')
                warning('''collchans'' option is deprecated'); 
            end          
            
            Sout.timewin = [-Inf Inf];
            Sout.phase   = S.tf.phase;
            Sout.method  = 'morlet';

            Sout.settings.subsample   = 1;
            Sout.settings.ncycles     = S.tf.Mfactor;
        else
            Sout = S;
        end
    otherwise
        Sout = S;
end
            
            
            
            
            