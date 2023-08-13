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
if ~isfield(MDP,'b'), MDP.b = MDP.B;  end
if ~isfield(MDP,'a'), MDP.a = MDP.A;  end


for f = 1:numel(MDP.b)

    % plot priors (transition probabilities)
    %----------------------------------------------------------------------
    for u = 1:size(MDP.b{f},3)
        subplot(6,6,u + (f - 1)*6)
        if strcmp(OPT,'norm')
            imagesc(spm_dir_norm(MDP.b{f}(:,:,u)))
        else
            imagesc(MDP.b{f}(:,:,u))
        end
        title('Transition priors')
        axis square
    end
end

% plot likelihood mapping
%--------------------------------------------------------------------------
subplot(2,2,3)
for g = 1:numel(MDP.a)
    a{g} = spm_dir_norm(MDP.a{g}(:,:));
end
A = spm_cat(a');
if size(A,1) > 256
    spm_spy(A);
else
    imagesc(A)
end
axis square, title('Likelihood','FontSize',14)
xlabel('latent states'),ylabel('outcomes')

% And allowable state transitions
%--------------------------------------------------------------------------
subplot(2,2,4)
imagesc(sum(MDP.b{1},3) > 1/32),axis square
title('Allowable transitions','FontSize',14), axis square
xlabel('latent states'),ylabel('latent states')