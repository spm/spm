function [R, RPV] = spm_mesh_resels(M,T,S,ndf)
% Returns the RESEL counts of a search volume on a surface mesh
% FORMAT R = spm_mesh_resels(M,T,[S])
% M        - a patch structure or [nx3] faces array (#faces = n)
% T        - a [mx1] logical vector (#vertices = m) defining search volume
% S        - a [mxp] array of standardised residuals [optional]
% ndf      - a 2-vector, [n df], the original n & dof of the linear model
%
% R        - a [1xD] array of RESEL counts {adimensional}
% RPV      - a [mx1] vector of RESELs per vertex
%__________________________________________________________________________
%
% References:
%
% [1] Detecting Sparse Signals in Random Fields, With an Application to
% Brain Imaging, J.E. Taylor and K.J. Worsley, Journal of the American
% Statistical Association, 102(479):913-928, 2007.
%
% [2] SurfStat: http://www.math.mcgill.ca/keith/surfstat/, K.J. Worsley.
%__________________________________________________________________________
% Copyright (C) 2010-2016 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_resels.m 6678 2016-01-14 18:23:33Z guillaume $

%-Parse input arguments
%--------------------------------------------------------------------------
if nargin < 3 || isempty(S)
    if ~isnumeric(M)
        S = M.vertices;
    else
        S = ones(max(M(:)),2);
    end
end
if nargin < 4
    ndf = [size(S,2) size(S,2)]; % Assume full df
end

if ~isnumeric(M), M = M.faces; end

if nargin < 2
    T = true(size(S,1),1);
end
if any(isnan(T)), T = ~isnan(T); end

D = 2;

%-Compute edges
%--------------------------------------------------------------------------
M      = sort(M,2); % so that ismember(M...,E,'rows') works later
E      = spm_mesh_edges(M);

%-Compute (standardised) residual sum of squares of differences along edges
%--------------------------------------------------------------------------
SSR    = S(E',:);
SSR    = reshape(SSR',size(SSR,2),2,[]);
SSR    = mean(squeeze((SSR(:,1,:) - SSR(:,2,:)).^2),1)';
SSR    = SSR * (ndf(1)/ndf(2)); % see comment in spm_est_smoothness.m
SSR    = sqrt(SSR);

%-Lipschitz-Killing Curvature (LKC)
%==========================================================================
LKC    = zeros(1,3);

%-Vertices contribution
%--------------------------------------------------------------------------
LKC(1) = LKC(1) + sum(T);

%-Edges contribution
%--------------------------------------------------------------------------
TE     = all(T(E),2);

LKC(1) = LKC(1) - sum(TE);
LKC(2) = LKC(2) + sum(SSR(TE));

%-Triangles contribution
%--------------------------------------------------------------------------
TM     = all(T(M),2);
% [E,I,J] = unique([M(:,[1 2]);M(:,[2 3]);M(:,[1 3])],'rows');
% SSR = reshape(SSR(J),[],3); SSR = SSR(TM,:); A = spm_mesh_area(SSR',true);
[is,l] = ismember(M(TM,[1 2]),E,'rows'); l12 = SSR(l);
[is,l] = ismember(M(TM,[1 3]),E,'rows'); l13 = SSR(l);
[is,l] = ismember(M(TM,[2 3]),E,'rows'); l23 = SSR(l);
A      = spm_mesh_area([l12 l13 l23]',true);

LKC(1) = LKC(1) + sum(TM);
LKC(2) = LKC(2) - sum(l12 + l13 + l23) / 2;
LKC(3) = LKC(3) + sum(A);

%-RESELS per vertex
%==========================================================================
if nargout > 1
    RPV = zeros(size(T));
    for i=1:3
        RPV = RPV + accumarray(M(TM,i),A,[size(T,1) 1]);
    end
    RPV = RPV / 3; % each triangle contributes a third to each of its vertices 
                   % (sum(A)==sum(RPV))
    RPV = RPV / sqrt((4*log(2))^D);
end

%-RESEL Counts
%==========================================================================
R      = LKC ./ sqrt((4*log(2)).^(0:D));
