function spm_MDP_params(MDP,OPT)
% Show likelihood and prior transition parameters
% FORMAT spm_MDP_params(MDP,OPT)
%
% MDP - MDP structure
% OPT - normalise [1]
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% deal with a sequence of trials
%==========================================================================
if nargin < 2,        OPT   = 'norm'; end
if ~isfield(MDP,'a'), MDP.a = MDP.A;  end
if ~isfield(MDP,'b'), MDP.b = MDP.B;  end


[Nf,Ns,Nu] = spm_MDP_size(MDP);

for f = 1:Nf

    % plot priors (transition probabilities)
    %----------------------------------------------------------------------
    for u = 1:Nu(f)
        subplot(6,6,u + (f - 1)*6)
        if strcmp(OPT,'norm')
            if Ns(f) > 128
                spm_spy(spm_dir_norm(MDP.b{f}(:,:,u)) > 1/16,4);
            else
                imagesc(spm_dir_norm(MDP.b{f}(:,:,u)))
            end
        else
            if Ns(f) > 128
                spm_spy(MDP.b{f}(:,:,u) > 1/16,4);
            else
                imagesc(MDP.b{f}(:,:,u))
            end
            imagesc(MDP.b{f}(:,:,u))
        end
        title('Transition priors')
        axis square
    end
end

% plot likelihood mapping
%--------------------------------------------------------------------------
try
    for g = 1:numel(MDP.a)
        a{g} = spm_dir_norm(MDP.a{g}(:,:));
    end
    A = spm_cat(a');

    subplot(2,2,3)
    if size(A,1) > 128
        spm_spy(A,8);
    else
        imagesc(A)
    end
    axis square, title('Likelihood','FontSize',14)
    xlabel('latent states'),ylabel('outcomes')
end


% And allowable state transitions
%--------------------------------------------------------------------------
subplot(2,2,4)
if Ns(f) > 128
    spm_spy(sum(MDP.b{1},3) > 1/16,8);
else
    imagesc(sum(MDP.b{1},3) > 1/16),
end
axis square
title('Allowable transitions','FontSize',14), axis square
xlabel('latent states'),ylabel('latent states')

