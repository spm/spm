function MDP = spm_set_costs(MDP,S,chi)
% Sets prior preferences and costs MDP{n}.C{g} in a renormalising MDP
% FORMAT MDP = spm_set_costs(MDP,S,chi)
% MDP  - Generative model (hierarchical)
% S    - list of modality streams; e.g. [2,3]
% chi  - log prior for each modality; e.g., [-8,8]  
%
%--------------------------------------------------------------------------
% This auxiliary routine generates outcomes at successive levels of a
% renormalising generative model to identify rewarding or costly outcomes
% at the lowest level, to set prior preferences or costs in MDP{n}.C{g}.
% These priors are specified in terms of chi, the log probability. This
% means a negative chi specifies a cost and a positive chi specifies a
% reward in the respective streams.
%
% Under these models,, a preferred or costly outcome corresponds to one or
% more outcomes generated within each stream. This means that preferences
% and costs are usually specified for distinct streams. In the example
% above, stream 2 reports costly outcomes, while stream 3 reports preferred
% outcomes.
%__________________________________________________________________________


% if multiple streams replace C (and U)
%==========================================================================
Nm    = numel(MDP);
Ns    = numel(S);
if Ns > 1
    for n = 1:Nm
        if isfield(MDP{n},'C')
            MDP{n} = rmfield(MDP{n},'C');
        end
    end
    for s = 1:Ns
        MDP = spm_set_costs(MDP,S(s),chi(s));
    end
    return

end

% Get sequences of generalised outcomes for stream S
%==========================================================================

% for each level, catagorise states of S in terms of outcomes
%--------------------------------------------------------------------------
C     = {};
for m = 1:(Nm - 1)
    
    % for each state at this level
    %----------------------------------------------------------------------
    sf     = MDP{m}.sB == S;              % factors of stream S 
    s      = [];                          % sequences of states
    c      = [];                          % contraints on outcomes
    for si = 1:size(MDP{m}.b{sf},1)

        % propagate predictions under this state down hierarchy
        %------------------------------------------------------------------
        s{m}  = si;
        for n = m:-1:2
            
            % accumate subordinate states of S
            %--------------------------------------------------------------
            s{n - 1} = [];
            sf    =  MDP{n - 1}.sB == S;
            sg    = find(MDP{n}.sA == S);
            for t = 1:size(s{n},2)

                % generate level n outcomes
                %----------------------------------------------------------
                [~,x] = max(MDP{n}.a{sg(1)}(:,s{n}(t)));
                [~,u] = max(MDP{n}.a{sg(2)}(:,s{n}(t)));

                % iterate over time under generated path
                %----------------------------------------------------------
                for r = 2:MDP{n}.T
                    [~,j] = max(MDP{n - 1}.b{sf}(:,x(r - 1),u));
                    x(r)  = j;
                end

                % accumulate path
                %----------------------------------------------------------
                s{n - 1} = [s{n - 1} x];
            end
        end

        % final outcomes
        %------------------------------------------------------------------
        o     = [];
        for t = 1:size(s{1},2)
            sg    = MDP{1}.sA == S;
            [~,j] = max(MDP{1}.a{sg}(:,s{1}(1,t)));
            o(t)  = j;
        end

        % does the outcome for state si contain a constraint?
        %------------------------------------------------------------------
        c(si,1) = any(o > 1);

    end

    % save contraints on the states of level m of stream S
    %----------------------------------------------------------------------
    C{m} = c;

end


% Place contraints in predictions (i.e., outcomes) of initial states
%==========================================================================
for n = 2:Nm

    % control factors and contraints
    %----------------------------------------------------------------------
    Ng    = numel(MDP{n}.a);                 % number of modalities
    Nf    = numel(MDP{n}.b);                 % number of factors

    if ~isfield(MDP{n},'U')
        MDP{n}.U = false(1,Nf);
    end
    if ~isfield(MDP{n},'C')
        for g = 1:Ng
            MDP{n}.C{g,1} = spm_dir_norm(ones(size(MDP{n}.a{g},1),1));
        end
    end

    % parents of stream S intial states
    %----------------------------------------------------------------------
    pf = MDP{n - 1}.id.D{MDP{n - 1}.sB == S};

    % states predicted by first stream
    %----------------------------------------------------------------------
    ps = find(ismember(MDP{n}.sA,1));

    % initial states of stream S predicted by first stream
    %----------------------------------------------------------------------
    pg = intersect(ps,pf);
    for g = pg
       
       % enable this factor
       %-------------------------------------------------------------------
       if spm_MDP_MI(MDP{n}.a{g}) > 1/512
           f = MDP{n}.id.A{g};
           MDP{n}.U(f) = true;
       end

       % and implement contraints on its predictions
       %-------------------------------------------------------------------
       MDP{n}.C{g} = spm_softmax(C{n - 1},chi);

    end

end

return
