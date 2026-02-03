function spm_RDP_params(MDP)
% Show prior transition parameters for a recursive or renormalising model
% FORMAT spm_RDP_params(MDP)
% MDP - MDP structure
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% deal with a sequence of trials
%==========================================================================

% unpack RDP if necessary
%--------------------------------------------------------------------------
if ~iscell(MDP)
    MDP = spm_rdp2mdp(MDP);
end

Nn    = numel(MDP);
for n = 1:Nn

    % plot priors (transition probabilities)
    %----------------------------------------------------------------------
    if ~isfield(MDP{n},'b')
        b = MDP{n}.B;
    else
        b = MDP{n}.b;
    end

    % plot priors (transition probabilities)
    %----------------------------------------------------------------------
    Nf    = numel(b);
    for f = 1:Nf
        B    = spm_dir_norm(sum(b{f},3));
        subplot(Nn,Nf,f + Nf*(Nn - n))
        if size(B,1) > 128
            spm_spy(B > 1/16,4);
        else
            imagesc(-B)
        end
        if Nf < 4, title('Transition priors'), end
        axis square
    end

end

