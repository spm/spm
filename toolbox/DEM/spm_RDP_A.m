function MDP = spm_RDP_A(MDP,PDP)
% places parameters in a recursive model
% FORMAT RDP = spm_RDP_A(MDP,PDP)
% MDP - recursive MDP
% PDP - recursive MDP
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

% Assume time scaling with a scale doubling
%--------------------------------------------------------------------------
MDP.a = PDP.a;
NL    = numel(PDP.Q.a);
for L = 1:NL
    str = repmat('MDP.',1,L);
    a   = PDP.Q.a{NL - L + 1};
    eval([str 'MDP.a = a;']) 
end

return