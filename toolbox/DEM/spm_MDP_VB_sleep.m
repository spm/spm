function [MDP] = spm_MDP_VB_sleep(MDP,OPTIONS)
% Bayean model erduction (sleep) for MDP models
% FORMAT [MDP] = spm_MDP_VB_sleep(MDP,OPTIONS)
%
% MDP(1:N) - cell array of (inverted) MDP stuctures
%
% See also: spm_MDP_VB and spm_MDP_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_VB_sleep.m 6747 2016-03-12 11:33:11Z karl $


% deal with a sequence of trials
%==========================================================================

% options
%--------------------------------------------------------------------------
try, OPTIONS.plot;    catch, OPTIONS.plot    = 0; end

% Baysian model reduction - A parameters
%--------------------------------------------------------------------------
g = 1;
a = MDP(1).a{g};
A = MDP(end).a{g};


% column-wise reduction
%--------------------------------------------------------------------------
sA     = A;
for i1 = 1:size(A,2)
    for i2 = 1:size(A,3)
        for i3 = 1:size(A,4)
            for i4 = 1:size(A,5)
                
                pA    = a(:,i1,i2,i3,i4);
                qA    = A(:,i1,i2,i3,i4);
                for i = 1:size(A,1);
                    rA    = pA;
                    rA(i) = 1;
                    dA    = qA + rA - pA;
                    dF    = spm_betaln(qA) + spm_betaln(rA) - spm_betaln(pA) - ...
                        spm_betaln(dA);
                    if dF < - 1/512;
                        sA(i,i1,i2,i3,i4) = 1;
                    end
                    
                end
                
            end
        end
    end
end


