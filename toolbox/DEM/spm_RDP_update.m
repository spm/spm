function MDP = spm_RDP_update(MDP,PDP,BMR)
% places parameters in a recursive model
% FORMAT RDP = spm_RDP_update(MDP,PDP,[BMR])
% MDP - nested MDP: prior
% PDP - nested MDP: posterior (with PDP.Q)
% BMR - Bayesian model reduction {'SIMPLE','MI'} [default: none]
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

% preliminaries
%--------------------------------------------------------------------------
if nargin < 3, BMR = []; end

% Dirichlet parameters to update
%--------------------------------------------------------------------------
param = {'a','b','c','d','e','h'};

% final level
%--------------------------------------------------------------------------
for p = 1:numel(param)
    if isfield(PDP.Q,param{p})

        % posteriors and priors
        %------------------------------------------------------------------
        qp = PDP.(param{p});
        pp = MDP.(param{p});

        % Bayesian model reduction
        %------------------------------------------------------------------
        if numel(BMR)
            for g = 1:numel(qp)
                qp{g} = spm_MDP_VB_prune(qp{g},pp{g},0,0,[],BMR);
            end
        end

        % update
        %------------------------------------------------------------------
        MDP.(param{p}) = qp;

    end
end

% subordinate levels
%--------------------------------------------------------------------------
NL    = numel(PDP.Q.a);
for L = 1:NL
    str   = repmat('MDP.',1,L);
    for p = 1:numel(param)

        if isfield(PDP.Q,param{p})

            % posteriors and priors
            %--------------------------------------------------------------
            qp = PDP.Q.(param{p}){NL - L + 1};
            pp = eval([str 'MDP.(param{p});']);

            % Bayesian model reduction
            %--------------------------------------------------------------
            if numel(BMR)
                for g = 1:numel(qp)
                    qp{g} = spm_MDP_VB_prune(qp{g},pp{g},0,0,[],BMR);
                end
            end

            % update
            %--------------------------------------------------------------
            eval([str 'MDP.(param{p}) = qp;']);

        end
    end
end


return