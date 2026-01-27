function pdp   = spm_VB_XXX(pdp)
% A fast form of MDP inversion
% FORMAT pdp   = spm_VB_XXX(pdp)
% pdp - likelihood (A) and transition (B) tensors for this sequence (O)
%
% This subroutine is a stripped down version of the MDP inversion schemes
% that exploits situations in which each likelihood mapping is a matrix;
% i.e., the latent states are conditionally independent because each
% likelihood has only one parent (and there are no priors).
%__________________________________________________________________________

% sizes
%--------------------------------------------------------------------------
Nf    = numel(pdp.B);
Ns    = ones(1,Nf);
Nu    = ones(1,Nf);
for f = 1:Nf
    Ns(f) = size(pdp.B{f},2);
    Nu(f) = size(pdp.B{f},3);
end

% each likelihood has only one parent
%--------------------------------------------------------------------------
iA    = full(spm_cat(pdp.id.A));

% states
%--------------------------------------------------------------------------
pdp.X = cell(Nf,1);
for t = 1:pdp.T
    for f = 1:Nf
        if Ns(f) > 1
            L = 0;
            for g = find(ismember(iA,f))
                L = spm_log(pdp.A{g}'*pdp.O{g,t}) + L;
            end
            pdp.X{f}(:,t) = spm_softmax(L);
        else
            pdp.X{f}(:,t) = 1;
        end
    end
end

% paths
%--------------------------------------------------------------------------
if pdp.T > 1
    pdp.P = cell(Nf,1);
    for t = 2:pdp.T
        for f = 1:Nf
            if Nu(f) > 1
                U     = zeros(Nu(f),1);
                for u = 1:Nu(f)
                    U(u) = pdp.X{f}(:,t)'*pdp.B{f}(:,:,u)*pdp.X{f}(:,t - 1);
                end
                pdp.P{f}(:,t - 1) = spm_softmax(spm_log(U));
            else
                pdp.P{f}(:,t - 1) = 1;
            end
        end
    end
else
    pdp.P = repmat({1},Nf,1);
end

return