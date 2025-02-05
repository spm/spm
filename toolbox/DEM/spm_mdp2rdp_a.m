function RDP = spm_mdp2rdp_a(MDP,p,q,T,FIX)
% Converts a cell array of MDPs into a recursive MDP (Dirichlet version)
% FORMAT RDP = spm_mdp2rdp_a(MDP,p,q,T,FIX)
% MDP{n} - Cell array of MDPs
% 
%  MDP{n}.a    - likelihood tensors (Dirichlet)
%  MDP{n}.b    - transition tensors (Dirichlet)
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

% Check for concentration parameters for streams
%--------------------------------------------------------------------------
Nm    = numel(MDP);
Ns    = max(MDP{1}.sB);
try p = p(1:Ns); catch, p = repmat(p(1),1,Ns); end
try q = q(1:Ns); catch, q = repmat(q(1),1,Ns); end

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
if numel(MDP{n}.b) > 1

    % remove trailing factors
    %----------------------------------------------------------------------
    MDP{n}.b = MDP{n}.b(1);
    MDP{n}.G = MDP{n}.G(1);

    if isfield(MDP{n},'sB')
        MDP{n}.sB = MDP{n}.sB(1);
    end

    Na    = numel(MDP{n}.a);
    d     = false(1,Na);
    for g = 1:Na
        d(g) = ~any(MDP{n}.id.A{g} > 1);
    end
    i = find(d);

    % remove their children
    %----------------------------------------------------------------------
    MDP{n}.a     = MDP{n}.a(d);
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
    d     = true(1,numel(MDP{n}.a));
    for g = 1:numel(MDP{n}.a)
        if (size(MDP{n}.a{g},1) < 2) && ~isa(MDP{n}.a{g},'function_handle')
            d(g) = false;
        end
    end
    i = find(d);

    % remove unitary likelihood mappings
    %----------------------------------------------------------------------
    MDP{n}.a    = MDP{n}.a(d);
    MDP{n}.id.A = MDP{n}.id.A(d);

    if isfield(MDP{n},'C')
        MDP{n}.C = MDP{n}.C(d);
    end

    if isfield(MDP{n},'sA')
        MDP{n}.sA = MDP{n}.sA(d);
        MDP{n}.sC = MDP{n}.sC(d);
    end

    % and update parents of subordinate factors
    %----------------------------------------------------------------------
    for j = 1:numel(MDP{n - 1}.id.D)
        MDP{n - 1}.id.D{j} = find(ismember(i,MDP{n - 1}.id.D{j}));
    end
    for j = 1:numel(MDP{n - 1}.id.E)
        MDP{n - 1}.id.E{j} = find(ismember(i,MDP{n - 1}.id.E{j}));
    end

    % remove unitary transitions
    %======================================================================
    d     = true(1,numel(MDP{n}.b));
    for f = 1:numel(MDP{n}.b)
        if isscalar(MDP{n}.b{f}) && ~isa(MDP{n}.b{f},'function_handle')
            d(f) = false;
        end
    end

    % leading constant factor
    %----------------------------------------------------------------------
    c    = find(~d,1,'first');
    d(c) = true;
    i    = find( d);
    k    = find(~d);

    % merge redundant factors (i.e., background states)
    %----------------------------------------------------------------------
    MDP{n}.b    = MDP{n}.b(d);
    MDP{n}.id.D = MDP{n}.id.D(d);
    MDP{n}.id.E = MDP{n}.id.E(d);

    if isfield(MDP{n},'U')
        MDP{n}.U = MDP{n}.U(d);
    end
    if isfield(MDP{n},'sB')
        MDP{n}.sB = MDP{n}.sB(d);
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
        Nf   = numel(MDP{n}.b);
        MDP{n}.U = false(1,Nf);
    end
end

% fill in empty columns of transition tensors (for controlled factors)
%==========================================================================
for n = 1:Nm
    for f = 1:numel(MDP{n}.b)

        % if this path is selectable (controllable)
        %------------------------------------------------------------------
        if MDP{n}.U(f)
            
            b  = MDP{n}.b{f};
            Ns = size(b,2);
            Nu = size(b,3);

            % preclude ambiguous transitions
            %--------------------------------------------------------------
            if Nu > 1 && Ns > 1 && ~isa(b,'function_handle')
                for u = 1:Nu
                    for s = 1:Ns
                        if ~any(b(:,s,u))
                            [j,i]    = max(max(squeeze(b(:,s,:)),[],2));
                            b(i,s,u) = j;
                        end
                    end
                end
            end
            MDP{n}.b{f} = b;

        end
    end
end


% prior concentration parameters
%==========================================================================
for n = 1:Nm

    % likelihoods
    %----------------------------------------------------------------------
    if FIX.A

        % normalise
        %------------------------------------------------------------------
        for g = 1:numel(MDP{n}.a)
            if ~isa(MDP{n}.a{g},'function_handle')
                MDP{n}.A{g} = spm_dir_norm(MDP{n}.a{g});
            end
        end
        MDP{n} = rmfield(MDP{n},'a');

    else

        % add concentration paramter
        %------------------------------------------------------------------
        for g = 1:numel(MDP{n}.a)
            if ~isa(MDP{n}.a{g},'function_handle')
                s           = MDP{n}.sC(g);
                MDP{n}.a{g} = MDP{n}.a{g} + p(s);
            end
        end

    end

    % priors
    %----------------------------------------------------------------------
    if FIX.B

        % normalise
        %------------------------------------------------------------------
        for f = 1:numel(MDP{n}.b)
            if ~isa(MDP{n}.b{f},'function_handle')
                s     = MDP{n}.sB(f);
                b     = spm_dir_norm(plus(MDP{n}.b{f},q(s)));

                MDP{n}.B{f} = b;
            end
        end
        MDP{n} = rmfield(MDP{n},'b');

    else

        % add concentration paramter
        %------------------------------------------------------------------
        for f = 1:numel(MDP{n}.b)
            if ~isa(MDP{n}.b{f},'function_handle')
                s           = MDP{n}.sB(f);
                MDP{n}.b{f} = MDP{n}.b{f} + q(s);
            end
        end
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


