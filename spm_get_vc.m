function SPM = spm_get_vc(SPM)
% generate variance components for factorial designs
% FORMAT SPM = spm_get_vc(SPM)
%
% SPM - SPM struct
% required fields
% SPM.xVi.I           - matrix containing levels for each scan and factor
% SPM.factor.variance - for each factor, indicate whether variances are
%                       equal (0) or unequal (1) between levels
% SPM.factor.dept     - for each factor, indicate whether variances are
%                       independent (0) or dependent (1) between levels
% set in output
% SPM.xVi.Vi          - cell vector of covariance components
%_______________________________________________________________________
%
% spm_get_vc generates variance components for a given design. For each
% factor, the user specifies whether its levels have identical variances
% and are uncorrelated. The individual components for each factor are
% combined into covariance components by using the Kronecker tensor
% product. If there are unequal number of observations at different levels,
% the function specifies covariance components for a full factorial
% design first and subsequently removes unwanted rows and columns from
% the covariance matrices.
%
% The functionality of spm_get_vc is similar to that of
% spm_non_sphericity. The difference is that spm_get_vc can accommodate 
% any number of factors and is more general, because it can cope with
% different number of observations under different levels of a factor.
%_______________________________________________________________________
% Copyright (C) 2006 Freiburg Brain Imaging 
% This code is part of SPM, which is
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Volkmar Glauche
% $Id: spm_get_vc.m 3468 2009-10-15 18:59:38Z karl $
 
% set up (numbers of scans and factors)
%--------------------------------------------------------------------------
Iin             = SPM.xVi.I;
[nscan nfactor] = size(Iin);
 
% make sure each row of Iin is unique
%==========================================================================
[Iu Ii Ij] = unique(Iin,'rows');
if size(Iu,1) < nscan
    nfactor = nfactor+1;
    uf = zeros(nscan, 1);
    for k = 1:max(Ij)
        uf(Ij==k) = 1:sum(Ij==k);
    end;
    Iin = [Iin uf];
end;
Nlevels = max(Iin);
Vi = {};
 
% first factor in SPM is replications, assume identical variance and independence
% pad with zeroes in case there are less than nfactor factors specified
%--------------------------------------------------------------------------
variance = [0 cat(2, SPM.factor.variance) zeros(1,nfactor)];
dept = [0 cat(2, SPM.factor.dept) zeros(1,nfactor)];
 
 
% (i) generate generic index
%==========================================================================
ngen = prod(Nlevels);
Igen = zeros(ngen, nfactor);
Igen(:,1) = kron(ones(1,prod(Nlevels(2:end))),1:Nlevels(1))';
for cf = 2:(nfactor-1)
    Igen(:,cf) = kron(ones(1,prod(Nlevels((cf+1):end))),kron(1:Nlevels(cf),ones(1,prod(Nlevels(1:(cf-1))))))';
end;        
Igen(:,nfactor) = kron(1:Nlevels(nfactor),ones(1,prod(Nlevels(1:(nfactor-1)))))';
        
% (ii) generate error variance components
%==========================================================================
for f=1:nfactor
    
    % identical/non-identical variances
    % for each factor, create a single variance component if variances are
    % identical across levels, and level specific variance components if
    % variances are non-identical
    %----------------------------------------------------------------------
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
    %----------------------------------------------------------------------
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
 
%(iii) sort out rows/columns & remove all-zero variance components
%==========================================================================
[unused ind] = ismember(Iin,Igen,'rows');
az = false(size(Vi));
 
for cVi = 1:numel(Vi)
    Vi{cVi} = Vi{cVi}(ind,ind);
    az(cVi) = full(all(Vi{cVi}(:) == 0));
end;
Vi = Vi(~az);
 
dupl = false(size(Vi));
for cVi = 1:numel(Vi)
    if ~dupl(cVi)
        for cVi1 = (cVi+1):numel(Vi)
            dupl(cVi1) = dupl(cVi1)||full(all(Vi{cVi}(:) == Vi{cVi1}(:)));
        end
    end;
end;
Vi = Vi(~dupl);
 
% (iv) save covariance components. If only one left, use this as error
% covariance matrix without going through the ReML step
%==========================================================================
if numel(Vi) == 1
    SPM.xVi.V = Vi{1};
else
    SPM.xVi.Vi = Vi;
end
