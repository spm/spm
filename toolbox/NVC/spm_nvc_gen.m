function [y] = spm_nvc_gen(P,M,U)
% Generate a BOLD signal prediction from scaled summed of neuronal drives 
% (neurovascular coupling).
% FORMAT [y] = spm_nvc_gen(P,M,U)
%
% Inputs:
% -------------------------------------------------------------------------
%  P - parameters of neurovascular coupling and Extended Balloon model
%  M - Neural mass model structure (M.input - neuronal drive functions)
%  U - Inputs
%
% Outputs:
% -------------------------------------------------------------------------
%  y - BOLD predictions
%
% This code scales neuronal drive signals by neurovascular coupling parameters
% and uses it as a single input (per each region) to a haemodynamic function.
% The outputs of the code are BOLD responses.
%__________________________________________________________________________
% Jafarian, A., Litvak, V., Cagnan, H., Friston, K.J. and Zeidman, P., 2019.
% Neurovascular coupling: insights from multi-modal dynamic causal modelling
% of fMRI and MEG. arXiv preprint arXiv:1903.07478.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Amirhossein Jafarian
% $Id: spm_nvc_gen.m 7713 2019-11-25 16:00:34Z spm $

% Neurovascular coupling signal (scaled summed neuronal drives)
%--------------------------------------------------------------------------
pf = M.input.input;
n = M.l;
a = [];

% Pre-synaptic
%--------------------------------------------------------------------------
if (strcmp(M.Model(1), 'pre'))
    if (M.input.num == 1  && size(M.pE.J,2) == n)
        for i = 1:size(pf,1)
            a(:,i) = squeeze(pf(i,:,:))'*P.J(:,i);
        end
    elseif (M.input.num == 1  && size(M.pE.J,2)== 1)
        for i = 1:size(pf,1)
            a(:,i) = squeeze(pf(i,:,:))'*P.J(1:end,1);
        end
    end
end

% Decomposed
%--------------------------------------------------------------------------
if (strcmp(M.Model(1), 'de'))
    if (M.input.num == 3  && size(M.pE.J,2)== 3 && numel(size(M.pE.J))==2)
        for i = 1:size(pf{1,1},1)
            a(:,i) = squeeze(pf{1,1}(i,:,:))'*P.J(:,1) + ...
                     squeeze(pf{1,2}(i,:,:))'*P.J(:,2) + ...
                     squeeze(pf{1,3}(i,:,:))'*P.J(:,3) ;
        end
        
    elseif (M.input.num == 2  && size(M.pE.J,2)== 2  && numel(size(M.pE.J))==2)
        for i = 1:size(pf{1,1},1)
            a(:,i) = squeeze(pf{1,1}(i,:,:))'*P.J(:,1) + ...
                     squeeze(pf{1,2}(i,:,:))'*P.J(:,2)  ;
        end
        
    elseif (M.input.num == 2  && size(M.pE.J,3)== n)
        for i = 1:size(pf{1,1},1)
            a(:,i) = squeeze(pf{1,1}(i,:,:))'*squeeze(P.J(:,1,i)) + ...
                     squeeze(pf{1,2}(i,:,:))'*squeeze(P.J(:,2,i)) ;
        end
    elseif (M.input.num == 3  && size(M.pE.J,3)== n)
        for i = 1:size(pf{1,1},1)
            a(:,i) = squeeze(pf{1,1}(i,:,:))'*squeeze(P.J(:,1,i)) + ...
                     squeeze(pf{1,2}(i,:,:))'*squeeze(P.J(:,2,i)) + ...
                     squeeze(pf{1,3}(i,:,:))'*squeeze(P.J(:,3,i)) ;
        end
    end
end

% Post-synaptic
%--------------------------------------------------------------------------
if (strcmp(M.Model(1), 'post'))
    if (M.input.num == 4  && size(M.pE.J,2)== n)
        for i = 1:size(pf,1)
            a(:,i) = squeeze(pf{i,1})'*P.J(1,i) + ...
                squeeze(pf{i,2})'*P.J(2,i) + ...
                squeeze(pf{i,3})'*P.J(3,i) + ...
                squeeze(pf{i,4})'*P.J(4,i) ;
        end
    elseif (M.input.num == 4  && size(M.pE.J,2)== 1)
        for i = 1:size(pf,1)
            a(:,i) = squeeze(pf{i,1})'*P.J(1,1) + ...
                squeeze(pf{i,2})'*P.J(2,1) + ...
                squeeze(pf{i,3})'*P.J(3,1) + ...
                squeeze(pf{i,4})'*P.J(4,1) ;
        end
    end
end

% Neurovascular signal
%--------------------------------------------------------------------------
U.u   = a;

% Haemodynamic model
%--------------------------------------------------------------------------
H.f   = @spm_fx_hdm;
H.g   = @spm_gx_hdm;
H.x   = M.x;
H.m   = M.m;
H.l   = M.l;
H.ns  = M.ns;

% Solve for haemodynamic responses
%--------------------------------------------------------------------------
y        =  spm_int(P.H,H,U);
