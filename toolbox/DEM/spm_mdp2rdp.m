function RDP = spm_mdp2rdp(MDP,p,q,T,FIX)
% Converts a cell array of MDPs into a recursive MDP
% FORMAT RDP = spm_mdp2rdp(MDP,p,q,T,FIX)
% MDP{n} - Cell array of MDPs
% 
%  MDP{n}.A    - likelihood tensors 
%  MDP{n}.B    - transition tensors (logical)
%  MDP{n}.id.A - cell array of parents of A factors at the same level
%  MDP{n}.id.D - cell array of parents of D in supraordinate outcomes
%  MDP{n}.id.E - cell array of parents of E in supraordinate outcomes
%
% p      - likelihood concentration; e.g., p = 1/32; [default: p = 0] 
% q      - prior (transition) decay; e.g., q = 1/32; [default: q = 0] 
% T      - path lengths [default: T = 2]
% 
% FIX.A = 0 for learning likelihoods
% FIX.B = 0 for learning transitions
%
% RDP    - likelihood and transition tensors for generating this sequence
%  RDP.L - level or depth
%  RDP.T - time steps
%
% This auxiliary routine takes a cell array of hierarchically arranged
% MDP’s and creates a single MDP where each subordinate MDP is a field of
% the superordinate MDP. The vertical dependencies are encoded in the cell
% arrays of parents at each level of the hierarchy, where the outcomes of
% one level are the parents of the initial states and paths of the lower
% level (encoded probabilistically in the vectors D and E of the lower
% level).
%
% In addition, this routine will remove redundant likelihood mappings, with
% only one output. Factors with only one state are consolidated into a
% single factor, which plays the role of a constant term; i.e., a unitary
% state that always generates the same outcome.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging



% Assume time scaling with a scale doubling
%--------------------------------------------------------------------------
if nargin < 2, p = 0; end
if nargin < 3, q = 0; end
if nargin < 4, T = 2; end
if nargin < 5
    FIX.A = 1;
    FIX.B = 1;
end

% is this RGM specified with Dirichlet counts?
%==========================================================================
if isfield(MDP{1},'a')
    RDP   = spm_mdp2rdp_a(MDP,p,q,T,FIX);
    return
end


% Check for concentration parameters
%--------------------------------------------------------------------------
Nm    = numel(MDP);
try p = p(1:Nm); catch, p = repmat(p(1),1,Nm); end
try q = q(1:Nm); catch, q = repmat(q(1),1,Nm); end


% Check for single level models
%--------------------------------------------------------------------------
if numel(MDP) < 2
    RDP   = MDP{1};
    RDP.L = 1;
    RDP.T = T;
    return
end

% if all streams comprise one group
%==========================================================================
n = Nm;
if numel(MDP{n}.B) > 1

    % remove trailing factors
    %----------------------------------------------------------------------
    MDP{n}.B = MDP{n}.B(1);
    MDP{n}.G = MDP{n}.G(1);

    Na    = numel(MDP{n}.A);
    d     = false(1,Na);
    for g = 1:Na
        d(g) = ~any(MDP{n}.id.A{g} > 1);
    end
    i = find(d);

    % remove their children
    %----------------------------------------------------------------------
    MDP{n}.A     = MDP{n}.A(d);
    MDP{n}.id.A  = MDP{n}.id.A(d);

    if isfield(MDP{n},'C')
        MDP{n}.C = MDP{n}.C(d);
    end

    % and update parents of subordinate factors
    %----------------------------------------------------------------------
    for j = 1:numel(MDP{n - 1}.id.D)
        MDP{n - 1}.id.D{j} = find(ismember(i,MDP{n - 1}.id.D{j}));
    end
    for j = 1:numel(MDP{n - 1}.id.E)
        MDP{n - 1}.id.E{j} = find(ismember(i,MDP{n - 1}.id.E{j}));
    end

end

% remove unitary mappings
%==========================================================================
for n = Nm:-1:2

    % find unitary likelihood mappings
    %======================================================================
    d     = true(1,numel(MDP{n}.A));
    for g = 1:numel(MDP{n}.A)
        if (size(MDP{n}.A{g},1) < 2) && ~isa(MDP{n}.A{g},'function_handle')
            d(g) = false;
        end
    end
    i = find(d);

    % remove unitary likelihood mappings
    %----------------------------------------------------------------------
    MDP{n}.A    = MDP{n}.A(d);
    MDP{n}.id.A = MDP{n}.id.A(d);

    if isfield(MDP{n},'C')
        MDP{n}.C = MDP{n}.C(d);
    end


    % and update parents of subordinate factors
    %----------------------------------------------------------------------
    for j = 1:numel(MDP{n - 1}.id.D)
        MDP{n - 1}.id.D{j} = find(ismember(i,MDP{n - 1}.id.D{j}));
    end
    for j = 1:numel(MDP{n - 1}.id.E)
        MDP{n - 1}.id.E{j} = find(ismember(i,MDP{n - 1}.id.E{j}));
    end

    % merge unitary transitions
    %======================================================================
    d     = true(1,numel(MDP{n}.B));
    for f = 1:numel(MDP{n}.B)
        if isscalar(MDP{n}.B{f}) && ~isa(MDP{n}.B{f},'function_handle')
            d(f) = false;
        end
    end

    % leading constant factor
    %----------------------------------------------------------------------
    c    = find(~d,1,'first');
    d(c) = true;
    i    = find( d);
    k    = find(~d);

    % remove redundant factors
    %----------------------------------------------------------------------
    MDP{n}.B    = MDP{n}.B(d);
    MDP{n}.id.D = MDP{n}.id.D(d);
    MDP{n}.id.E = MDP{n}.id.E(d);
    
    if isfield(MDP{n},'U')
        MDP{n}.U = MDP{n}.U(d);
    end

    % and update parents of likelihoods
    %----------------------------------------------------------------------
    for j = 1:numel(MDP{n}.id.A)
        if ismember(MDP{n}.id.A{j},k)
            MDP{n}.id.A{j} = c;
        else
            MDP{n}.id.A{j} = find(ismember(i,MDP{n}.id.A{j}));
        end
    end

end

% control factors
%-------------------------------------------------------------------------
for n = 1:Nm
    if ~isfield(MDP{n},'U')
        Nf   = numel(MDP{n}.B);
        MDP{n}.U = false(1,Nf);
    end
end

% fill in empty columns of transition tensors (for controlled factors)
%==========================================================================
for n = 1:Nm
    for f = 1:numel(MDP{n}.B)

        % if this path is selectable (controllable)
        %------------------------------------------------------------------
        if MDP{n}.U(f)
            
            b  = MDP{n}.B{f};
            Ns = size(b,2);
            Nu = size(b,3);

            % preclude ambiguous transitions
            %--------------------------------------------------------------
            if Nu > 1 && Ns > 1 && ~isa(b,'function_handle')
                for u = 1:Nu
                    for s = 1:Ns
                        if ~any(b(:,s,u))
                            i = find(any(squeeze(b(:,s,:)),2),1);
                            b(i,s,u) = true;
                        end
                    end
                end
            end
            MDP{n}.B{f} = b;

        end
    end
end


% prior concentration parameters
%==========================================================================
for n = 1:Nm

    % normalised likelihoods
    %----------------------------------------------------------------------
    for g = 1:numel(MDP{n}.A)
        if isa(MDP{n}.A{g},'function_handle')
            MDP{n}.a{g} = MDP{n}.A{g};
        else
            MDP{n}.a{g} = MDP{n}.A{g} + p(n);
        end
    end

    % normalised transitions
    %----------------------------------------------------------------------
    for f = 1:numel(MDP{n}.B)
        if isa(MDP{n}.B{f},'function_handle')
            MDP{n}.b{f} = MDP{n}.B{f};
        else
            MDP{n}.b{f} = MDP{n}.B{f} + q(n);
        end
    end

    % remove fields
    %----------------------------------------------------------------------
    if FIX.A
        MDP{n}.A = spm_dir_norm(MDP{n}.a);
        MDP{n}   = rmfield(MDP{n},'a');
    else
        MDP{n}   = rmfield(MDP{n},'A');
    end
    if FIX.B
        MDP{n}.B = spm_dir_norm(MDP{n}.b);
        MDP{n}   = rmfield(MDP{n},'b');
    else
        MDP{n}   = rmfield(MDP{n},'B');
    end

end


% remove parents of last level (i.e., MDP{end}.D and E)
%--------------------------------------------------------------------------
try
    MDP{end}.id = rmfield(MDP{end}.id,'D');
    MDP{end}.id = rmfield(MDP{end}.id,'E');
end

% Recursively place MDP in MDP.MDP, and specify path lengths (T)
%--------------------------------------------------------------------------
SDP   = MDP{1};
if ~isfield(SDP,'T')
    SDP.T = T;
end
SDP.L = 1;
for n = 2:Nm
    RDP     = MDP{n};
    RDP.MDP = SDP;
    SDP     = RDP;
    SDP.T   = T;
    SDP.L   = n;
end
RDP.L = n;

return


