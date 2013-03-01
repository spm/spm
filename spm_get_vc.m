function Vi = spm_get_vc(I,factor)
% Generate error covariance components for factorial designs
% FORMAT Vi = spm_get_vc(I,factor)
% I         - n x m matrix of factor level indicators
%             I(n,i) is the level of factor i for observation n
% factor(i) - structure array of sphericity assumptions for each factor
% .variance - 1 for different variance among levels of factor i
% .dept     - 1 for dependencies within levels of factor i
%
% Vi        - cell vector of covariance components
%__________________________________________________________________________
%
% spm_get_vc generates variance components for a given design. For each
% factor, the user specifies whether its levels have identical variances
% and are independent. The individual components for each factor are
% combined into covariance components by using the Kronecker tensor
% product. If there are unequal number of observations at different levels,
% the function specifies covariance components for a full factorial design
% first and subsequently removes unwanted rows and columns from the
% covariance matrices.
%
% The functionality of spm_get_vc is similar to that of spm_non_sphericity.
% The difference is that spm_get_vc can accommodate any number of factors
% and is more general, because it can cope with different number of
% observations under different levels of a factor.
%__________________________________________________________________________
% Copyright (C) 2006 Freiburg Brain Imaging 
% Copyright (C) 2008-2013 Wellcome Trust Centre for Neuroimaging
 
% Volkmar Glauche
% $Id: spm_get_vc.m 5293 2013-03-01 16:41:46Z guillaume $
 

%-Numbers of scans and factors
%--------------------------------------------------------------------------
[nscan,nfactor] = size(I);
 
%-Make sure each row of Iin is unique
%==========================================================================
[Iu,Ii,Ij]  = unique(I,'rows');
if size(Iu,1) < nscan
    nfactor = nfactor + 1;
    uf      = zeros(nscan, 1);
    for k = 1:max(Ij)
        uf(Ij==k) = 1:sum(Ij==k);
    end
    I       = [I uf];
end
nlevel      = max(I);
 
%-Non-sphericity assumptions
%--------------------------------------------------------------------------
% First factor is replications, assume identical variance and independence.
% Pad with zeroes in case there are less than nfactor factors specified.
variance    = [0 cat(2, factor.variance) zeros(1,nfactor)];
dept        = [0 cat(2, factor.dept) zeros(1,nfactor)];

% (i) generate generic index
%==========================================================================
Igen = zeros(prod(nlevel), nfactor);
Igen(:,1) = kron(ones(1,prod(nlevel(2:end))),1:nlevel(1))';
for cf = 2:(nfactor-1)
    Igen(:,cf) = kron(ones(1,prod(nlevel((cf+1):end))),kron(1:nlevel(cf),ones(1,prod(nlevel(1:(cf-1))))))';
end        
Igen(:,nfactor) = kron(1:nlevel(nfactor),ones(1,prod(nlevel(1:(nfactor-1)))))';
        
% (ii) generate error variance components
%==========================================================================
Vi = {};
for f=1:nfactor
    
    % identical/non-identical variances
    % for each factor, create a single variance component if variances are
    % identical across levels, and level specific variance components if
    % variances are non-identical
    %----------------------------------------------------------------------
    nVi = {};
    if ~variance(f)
        nVi{1} = speye(nlevel(f),nlevel(f));
    else
        for l1=1:nlevel(f)
            nVi{l1} = sparse(l1,l1,1,nlevel(f),nlevel(f));
        end
    end
    if dept(f)
        for l1 = 1:nlevel(f)
            for l2 = 1:(l1-1)
                nVi{end+1} = sparse([l1 l2],[l2 l1],1,nlevel(f), ...
                                         nlevel(f));
            end
        end
    end
    
    % combine current factor components with previous ones, thus building
    % up covariance components block by block
    %----------------------------------------------------------------------
    if isempty(Vi)
        Vi = nVi;
    else
        oVi = Vi;
        Vi = {};
        for nv = 1:numel(nVi)
            for ov = 1:numel(oVi)
                Vi{end+1} = kron(nVi{nv}, oVi{ov});
            end
        end
    end
end
 
% (iii) sort out rows/columns & remove all-zero variance components
%==========================================================================
[unused,ind] = ismember(I,Igen,'rows');
az = false(size(Vi));
 
for cVi = 1:numel(Vi)
    Vi{cVi} = Vi{cVi}(ind,ind);
    az(cVi) = full(all(Vi{cVi}(:) == 0));
end
Vi = Vi(~az);
 
dupl = false(size(Vi));
for cVi = 1:numel(Vi)
    if ~dupl(cVi)
        for cVi1 = (cVi+1):numel(Vi)
            dupl(cVi1) = dupl(cVi1)||full(all(Vi{cVi}(:) == Vi{cVi1}(:)));
        end
    end
end
Vi = Vi(~dupl);
