function con = spm_design_contrasts(SPM)
% Make contrasts for one, two or three-way ANOVAs
% FORMAT con = spm_design_contrasts(SPM)
% SPM           - SPM structure
%
% con           - structure array of contrasts with fields
%   con(i).c    - Contrast matrix
%   con(i).name - Name
%__________________________________________________________________________
%
% This function generates contrasts on the basis of the current SPM
% design. This is specified in SPM.factor (how the factors relate to the
% conditions) and SPM.xBF.order (how many basis functions per condition).
%
% This function generates (transposed) contrast matrices to test
% for the average effect of condition, main effects of factors and
% interactions.
%__________________________________________________________________________
% Copyright (C) 2005-2016 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_design_contrasts.m 6892 2016-09-30 12:56:35Z guillaume $


% Only relevant if factorial design has been specified
%--------------------------------------------------------------------------
if isempty(SPM.factor), con = []; return; end

% Compute contrasts for a generic k(1)-by-k(2)-by-k(3) design
%--------------------------------------------------------------------------
nf   = length(SPM.factor);
kf   = [SPM.factor.levels];

icon = spm_make_contrasts(kf);
ncon = length(icon);
if ~ncon, con = []; return; end

% Get number of basis functions per condition and 'expand' contrasts
%--------------------------------------------------------------------------
try
    nbases   = SPM.xBF.order;
catch
    nbases   = 1; % for 2nd level designs
end
for c=1:ncon
    con(c).c = kron(icon(c).c,eye(nbases));
end

% Handle parametric modulations/multi-sessions for 1st level designs
%--------------------------------------------------------------------------
if isfield(SPM,'Sess')

    nsess = length(SPM.Sess);
    conds = length(SPM.Sess(1).U);
    
    % Pad out parametric modulation columns with zeros
    %----------------------------------------------------------------------
    for c=1:ncon
        nr    = size(con(c).c,1);
        col   = 1;
        block = [];
        for cc=1:conds
            block = [block,con(c).c(:,col:col+nbases-1)];
            block = [block,zeros(nr,sum([SPM.Sess(1).U(cc).P.h]))];
            col   = col + nbases;
        end
        con(c).c  = block;
    end

    % If there are multiple sessions, replicate each contrast matrix over
    % sessions - this assumes the conditions are identical in each session
    %----------------------------------------------------------------------
    if nsess > 1
        covs  = zeros(1,nsess);
        for s=1:nsess
            nconds   = length(SPM.Sess(s).U);
            if nconds ~= conds
                error('Number of conditions must be same in all sessions');
            end
            covs(s)  = size(SPM.Sess(s).C.C,2);
        end
        for c=1:ncon
            nr       = size(con(c).c,1);
            c1       = [con(c).c zeros(nr,covs(1))];
            for s=2:nsess
                c1   = [c1 con(c).c zeros(nr,covs(s))];
            end
            con(c).c = c1;
        end
    end

    % Pad contrasts out - add columns for session effects
    %----------------------------------------------------------------------
    k = size(SPM.xX.X,2);
    for c=1:ncon
        con(c).c(:,end+1:k) = 0;
    end

end

% Rename contrasts (see spm_make_contrasts)
%--------------------------------------------------------------------------
con(1).name = 'Average effect of condition';
con(2).name = ['Main effect of ',SPM.factor(1).name];
if nf>1
    con(3).name = ['Main effect of ',SPM.factor(2).name];
    con(4).name = ['Interaction: ',SPM.factor(1).name,' x ',SPM.factor(2).name];
end
if nf>2
    con(5).name = ['Main effect of ',SPM.factor(3).name];
    con(6).name = ['Interaction: ',SPM.factor(1).name,' x ',SPM.factor(3).name];
    con(7).name = ['Interaction: ',SPM.factor(2).name,' x ',SPM.factor(3).name];
    con(8).name = ['Interaction: ',SPM.factor(1).name,' x ',SPM.factor(2).name,' x ',SPM.factor(3).name];
end
