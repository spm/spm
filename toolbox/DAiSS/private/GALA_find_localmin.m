function regions = GALA_find_localmin(lJcov,Nl,Nd,A,thresh)    

  
    for p=1:Nl           % loop by persons
    
        lJcovp{p}=lJcov(1+(p-1)*Nd:Nd+(p-1)*Nd);   % separate for individual person

        [slJcovp, islJcovp] = sort(lJcovp{p},'descend');
        Aip = A(islJcovp,islJcovp);             % reindex adjacency matrix
        slJcovp(slJcovp<=thresh*max(slJcovp)) = 0;               % set zero (and negative) values to -1
        
        % find boundary of basins for every local maximum
        
        nclb = 0;       % initialize number of blobs
        bound = [];
        
        Ndd = length(find(slJcovp>0));  % number of nonzero vertices
        
        for j=1:Ndd                     % loop by nonzero vertices
            
            % find vertices above slJcovp(j)
            up = spdiags(slJcovp'>=slJcovp(j),0,length(Aip),length(Aip));
            
            reg = up*Aip*up;    % adjacency matrix above slJcovp(j)
            [v_cl,cl_sz]=get_components(reg);   % get connected components
            
            cll=find(cl_sz >= 2);   % find indices of blobs (blob is >=2 anjacent vertices)
            
            % check if current vertex joins blobs
            % i.e. reduces number of blobs and is not disconnected (for
            % equal values)
            if length(cll)<nclb && cl_sz(v_cl(j))>1
                % then this vertex is boundary vertex
                bound = [bound; j];
                slJcovp(j) = 0; % set boundary vertex to zero
            else
                nclb=length(cll);   % else update number of blobs
                % we must not change the number of blobs for boundary vertex
                % because it removed (set to zero) from vertices above the
                % threshold and so nclb stay the same
            end          
        end
              
        % find indices of clusters and local maxima

        maxsp = [];

        up = spdiags(slJcovp'>0,0,length(Aip),length(Aip));
        reg = up*Aip*up;
        [v_cl,cl_sz]=get_components(reg);
        cllp{p}=find(cl_sz >= 2);
        
        clviA=[];
        for i=1:length(cllp{p})
            clvi{p,i} = find(v_cl==cllp{p}(i));
            [nu mi] = max(slJcovp(clvi{p,i}));
            maxsp = [maxsp; clvi{p,i}(mi)];
            clviA = [clviA clvi{p,i}];
        end 
        
        maxs{p} = maxsp;

        % add boundary vertices to appropriate cluster
        
        A0=Aip-eye(length(Aip));    % direct neighbors matrix
        for i=1:length(bound)   % loop by all boudary vertices
            
            neib = find(A0(:,bound(i)));    % direct neighbors of current verex

            cln = Inf;  % initialize number of appropriate cluster
            
            for j=1:length(neib)    % loop by all neighbors
                % find cluster index for current neighbor
                ind = find(cllp{p}==v_cl(neib(j)));
                % if there is no cluster set ind to Inf
                if isempty(ind);
                    ind = Inf;
                end
                % choose minimum of cluster numbers providing clusters with
                % bigger local maximum (less index) have priority
                cln = min(cln,ind);
            end
            % add vertex to the cluster with cln number
            clvi{p,cln}=sort([clvi{p,cln} bound(i)]);
        end
        clviAb = sort([clviA bound']);
        
        % reindex back all indices
        for i=1:length(cllp{p})
            regions.clvi{p,i} = islJcovp(clvi{p,i});
        end
        regions.boundp{p} = islJcovp(bound);
        regions.maxs{p} = islJcovp(maxs{p});
        regions.clviAp{p} = islJcovp(clviAb);

    end
