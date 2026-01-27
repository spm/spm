function MDP = spm_RDP_reduce(MDP,BMR,varargin)
% places parameters in a recursive model
% FORMAT RDP = spm_RDP_reduce(MDP,BMR)
% MDP   - nested MDP: full
% BMR   - Bayesian model reduction {'SIMPLE','MI','HARD','SOFT'}
% param - {'a','A',...}
%
% This routine applies various forms of Bayesian model reduction to the
% Dirichlet parameters or mappings specified in param (as strings). The
% simple or mutual information options will apply Bayesian model reduction
% based upon comparing models with and without a particular parameter
% (simple) or based upon models that increase mutual information (i.e.,
% decrease expected free energy). The hard and soft options simply
% threshold low probability entries; either using a hard threshold or by
% increasing the precision of the specified probabilistic mapping. Note
% that a soft option will convert Dirichlet parameters into a multinomial
% distribution (that may require rescaling if necessary).
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging


% parameters to update
%--------------------------------------------------------------------------
param = varargin;
beta  = 512;

% final level
%--------------------------------------------------------------------------
for p = 1:numel(param)
    if isfield(MDP,param{p})

        % posteriors and priors
        %------------------------------------------------------------------
        qp = MDP.(param{p});

        % Bayesian model reduction
        %------------------------------------------------------------------
        switch BMR

            case({'SIMPLE','MI'})
                for g = 1:numel(qp)
                    pp    = ones(size(qp{g}));
                    qp{g} = spm_MDP_VB_prune(qp{g},pp,0,0,0,BMR);
                end

            case({'HARD'})
                for g = 1:numel(qp)
                    i        = qp{g} < max(qp{g}(:))/16;
                    qp{g}(i) = 0;
                end

            case({'SOFT'})
                for g = 1:numel(qp)
                    qp{g} = spm_softmax(beta*spm_log(qp{g}));
                end

            otherwise
        end
        % update
        %------------------------------------------------------------------
        MDP.(param{p}) = qp;

    end
end

% subordinate levels
%--------------------------------------------------------------------------
NL    = MDP.L;
for L = 1:NL
    str   = repmat('MDP.',1,L);
    for p = 1:numel(param)

        try

            % posteriors and priors
            %--------------------------------------------------------------
            eval(['qp = ' str 'MDP.(param{p});']);

            % Bayesian model reduction
            %--------------------------------------------------------------
            switch BMR

                case({'SIMPLE','MI'})
                    for g = 1:numel(qp)
                        pp    = ones(size(qp{g}));
                        qp{g} = spm_MDP_VB_prune(qp{g},pp,0,0,0,BMR);
                    end

                case({'HARD'})
                    for g = 1:numel(qp)
                        i        = qp{g} < max(qp{g}(:))/16;
                        qp{g}(i) = 0;
                    end

                case({'SOFT'})
                    for g = 1:numel(qp)
                        qp{g} = spm_softmax(beta*spm_log(qp{g}));
                    end

                otherwise
            end

            % update
            %--------------------------------------------------------------
            eval([str 'MDP.(param{p}) = qp;']) 

        end
    end
end


return