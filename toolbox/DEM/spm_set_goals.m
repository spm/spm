function MDP = spm_set_goals(MDP,S,chi)
% Gets rewarded (and restricted) states at the deepest level of an MDP
% FORMAT MDP = spm_set_goals(MDP,S,chi)
% MDP  - Generative model (hierarchical)
% S    - list of modality streams;    e.g., [ 2,3]
% chi  - log prior for each modality; e.g., [-8,8]  
%
%--------------------------------------------------------------------------
% This auxiliary routine identifies the episodes (i.e., paths) encoded by
% states at the deepest or highest level of a hierarchical MDP that entail
% intended or avoided outcomes at the lowest level.These intended or costly
% latent states are indexed in MDP{end}.id.hid and MDP{end}.id.cid,
% respectively — in the highest level identifier field.
%__________________________________________________________________________


% if multiple streams replace C (and U)
%==========================================================================
Nm    = numel(MDP);
Ns    = numel(S);
if Ns > 1
    for s = 1:Ns
        MDP = spm_set_goals(MDP,S(s),chi(s));
    end
    return
end

% Get sequences of generalised outcomes for stream S
%==========================================================================
m = Nm;
if chi >= 0
    MDP{m}.id.hid = [];
    if isfield(MDP{m},'H')
        MDP{m}    = rmfield(MDP{m},'H');
    end
end
if chi <= 0
    MDP{m}.id.cid = [];
end

% parents of factor of stream S
%--------------------------------------------------------------------------
pd = MDP{m - 1}.id.D{MDP{m - 1}.sB == S};
pe = MDP{m - 1}.id.E{MDP{m - 1}.sB == S};

% inital states predicted by first factor of first stream
%--------------------------------------------------------------------------
i  = find(MDP{m}.sB == 1,1,'first');
ps = find(ismember([MDP{m}.id.A{:}],i));

% initial states of stream S predicted by first stream
%--------------------------------------------------------------------------
pd = intersect(ps,pd);
pe = intersect(ps,pe);

% for each state at this level
%--------------------------------------------------------------------------
s  = {};                                              % sequences of states
Ns = size(MDP{m}.b{1},1);
for si = 1:Ns

    % generate level m outcomes under this state
    %----------------------------------------------------------------------
    if numel(pd)
        [~,x] = max(MDP{m}.a{pd}(:,si));
        [~,u] = max(MDP{m}.a{pe}(:,si));
    else
        x = 1;
        u = 1;
    end

    % iterate over time under generated path
    %----------------------------------------------------------------------
    sf    =  MDP{m - 1}.sB == S;
    for r = 2:MDP{m}.T
        [~,j] = max(MDP{m - 1}.b{sf}(:,x(r - 1),u));
        x(r)  = j;
    end
    s{m - 1}  = x;

    % propagate predictions under this state down hierarchy
    %----------------------------------------------------------------------
    for n = flip(2:(m - 1))

        % accumate subordinate states of S
        %------------------------------------------------------------------
        s{n - 1} = [];
        sf    =  MDP{n - 1}.sB == S;
        sg    = find(MDP{n}.sA == S);
        for t = 1:size(s{n},2)

            % generate level n outcomes
            %--------------------------------------------------------------
            [~,x] = max(MDP{n}.a{sg(1)}(:,s{n}(t)));
            [~,u] = max(MDP{n}.a{sg(2)}(:,s{n}(t)));

            % iterate over time under generated path
            %--------------------------------------------------------------
            for r = 2:MDP{n}.T
                [~,j] = max(MDP{n - 1}.b{sf}(:,x(r - 1),u));
                x(r)  = j;
            end

            % accumulate path
            %--------------------------------------------------------------
            s{n - 1} = [s{n - 1} x];
        end
    end

    % final outcomes
    %----------------------------------------------------------------------
    o     = [];
    for t = 1:size(s{1},2)
        sg    = MDP{1}.sA == S;
        [~,j] = max(MDP{1}.a{sg}(:,s{1}(1,t)));
        o(t)  = j;
    end

    % Record hidden states generating this outcome
    %----------------------------------------------------------------------
    if any(o > 1)
        if chi > 0
            MDP{m}.id.hid(end + 1) = si;
        elseif chi < 0
            MDP{m}.id.cid(end + 1) = si;
        end
    end

end

% specify H
%--------------------------------------------------------------------------
if chi > 0 && numel(MDP{m}.id.hid)
    h = sparse(MDP{m}.id.hid,1,chi,Ns,1);
    MDP{m}.H{1} = spm_softmax(h);
end

% turn on control
%--------------------------------------------------------------------------
if chi > 0 && numel(MDP{m}.id.hid)
    MDP{m}.U = 1;
end
if chi < 0 && numel(MDP{m}.id.cid)
    MDP{m}.U = 1;
end

return
