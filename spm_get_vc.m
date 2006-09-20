function SPM=spm_get_vc(SPM)

Nlevels = max(SPM.xVi.I);
[nscan nfactor] = size(SPM.xVi.I);
Vi = {};

% first factor in SPM is replications, assume identical variance and independence
variance = [0 cat(2, SPM.factor.variance) zeros(1,nfactor)];
dept = [0 cat(2, SPM.factor.dept) zeros(1,nfactor)];

% first, generate generic index
Igen = zeros(prod(Nlevels), nfactor);
Igen(:,1) = kron(ones(1,prod(Nlevels(2:end))),1:Nlevels(1))';
for cf=2:(nfactor-1)
    Igen(:,cf) = kron(ones(1,prod(Nlevels((cf+1):end))),kron(1:Nlevels(cf),ones(1,prod(Nlevels(1:(cf-1))))))';
end;        
Igen(:,nfactor) = kron(1:Nlevels(nfactor),ones(1,prod(Nlevels(1:(nfactor-1)))))';
        
ngen = prod(Nlevels);

% second, generate error variance components
for f=1:nfactor
    % identical/non-identical variances
    % for each factor, create a single variance component if variances are
    % identical across levels, and level specific variance components if
    % variances are non-identical
    nVi = {};
    if ~variance(f)
        nVi{1} = speye(Nlevels(f),Nlevels(f));
    else
        for l1=1:Nlevels(f)
            nVi{l1} = sparse(l1,l1,1,Nlevels(f),Nlevels(f));
        end;
    end;
    if dept(f)
        for l1 = 1:Nlevels(f)
            for l2 = 1:(l1-1)
                nVi{end+1} = sparse([l1 l2],[l2 l1],1,Nlevels(f), ...
                                         Nlevels(f));
            end;
        end;
    end;
    % combine current factor components with previous ones, thus building
    % up covariance components block by block
    if isempty(Vi)
        Vi = nVi;
    else
        oVi = Vi;
        Vi = {};
        for nv = 1:numel(nVi)
            for ov = 1:numel(oVi)
                Vi{end+1} = kron(nVi{nv}, oVi{ov});
            end;
        end;
    end;
end;

% third, sort out rows/columns for real design
[unused ind] = ismember(SPM.xVi.I,Igen,'rows');

for cVi = 1:numel(Vi)
    Vi{cVi} = Vi{cVi}(ind,ind);
end;

% last, remove dupl & all-zero variance components
dupl = zeros(size(Vi));
for cVi = 1:numel(Vi)
    dupl(cVi) = dupl(cVi)||full(all(Vi{cVi}(:) == 0));
    if ~dupl(cVi)
        for cVi1 = (cVi+1):numel(Vi)
            dupl(cVi1) = dupl(cVi1)||full(all(Vi{cVi}(:) == Vi{cVi1}(:)));
        end
    end;
end;
Vi = Vi(~dupl);

% save covariance components. If only one left, use this as error
% covariance matrix without going through the ReML step
if numel(Vi) == 1
    SPM.xVi.V = Vi{1};
else
    SPM.xVi.Vi = Vi;
end;
