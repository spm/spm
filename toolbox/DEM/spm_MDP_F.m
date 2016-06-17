function [F,Fu,Fs,Fq,Fg,Fa] = spm_MDP_F(MDP)
% auxiliary function for retrieving free energy and its components
% FORMAT [F,Fu,Fs,Fq,Fg,Fa] = spm_MDP_F(MDP)
%
% F   - total free energy
% Fu  - confidence
% Fs  - free energy of states
% Fq  - free energy of policies
% Fg  - free energy of precision
% Fa  - free energy of parameters
%
% If MDP is a cell array, the free actions are turned (summed over time),
% otherwise, the free energies are turned over time
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_F.m 6811 2016-06-17 09:55:47Z karl $


% evaluate free action
%==========================================================================
m = numel(MDP);
if m > 1
    for i = 1:numel(MDP)
        
        % free action due to states and policies
        %------------------------------------------------------------------
        [f,fu,fs,fq,fg] = spm_MDP_F(MDP(i));
        F(i ) = sum(f);
        Fu(i) = sum(fu);
        Fs(i) = sum(fs);
        Fq(i) = sum(fq);
        Fg(i) = sum(fg);
        
        % free energy due to parameters
        %------------------------------------------------------------------
        try
            Fa = spm_cat({MDP.Fa});
        catch
            Fa = 0;
        end
    end
    return
    
else
    
    % evaluate free energy
    %======================================================================
    pg  = 1;                                    % prior precision
    qg  = MDP.w;                                % posterior precision
    pu  = spm_softmax(MDP.G*diag(qg));          % prior policies
    qu  = spm_softmax(MDP.F + MDP.G*diag(qg));  % posterior policies
    
    Fu  =  sum(qu.*log(qu));                    % confidence
    Fs  = -sum(qu.*MDP.F);                      % free energy of states
    Fq  = -sum(qu.*log(pu));                    % free energy of policies
    Fg  = qg/pg - log(qg);                      % free energy of precision
    Fa  = [];
    F   = Fs + Fu + Fq + Fg;                    % total free energy
    
end