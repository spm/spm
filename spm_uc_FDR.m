function [u,Ps,Ts] = spm_uc_FDR(q,df,STAT,n,Vs,Vm)
% False Discovery critical height threshold
% FORMAT [u,Ps,Ts] = spm_uc_FDR(q,df,STAT,n,Vs[,Vm])
%
% q     - critical expected False Discovery Rate
% df    - [df{interest} df{residuals}]
% STAT  - Statistical field  (see comments below about FWER and EFDR)
%       'Z' - Gaussian field
%       'T' - T - field
%       'X' - Chi squared field
%       'F' - F - field
%       'P' - P - value
% n     - Conjunction number
% Vs    - Mapped statistic image(s)
%          -or-
%         Vector of sorted p-values, p1<p2<... (saves i/o w/ repeated calls)
% Vm    - Mask in 1 of 3 forms
%           o Scalar, indicating implicit mask value in statistic image(s)
%           o Vector of indicies of elements within mask
%           o Mapped mask image
%         (Ignored if Vs is a vector.)
%
% u     - critical height
% Ps    - Sorted p-values
% Ts    - Sorted statistic values
%
%__________________________________________________________________________
%
% The Benjamini & Hochberch (1995) False Discovery Rate (FDR) procedure
% finds a threshold u such that the expected FDR is at most q.
% spm_uc_FDR returns this critical threshold u.
%
% For repeated use of a given statistic image, return Ps in the place
% of Vs:
%      [P Ps] = spm_uc_FDR(Z1,df,STAT,n,Vs);  %-Initialization, image read
%      P      = spm_uc_FDR(Z2,df,STAT,n,Ps);  %-Image not read, Ps used
%
% Note that a threshold of Inf is possible if q is very small.  This
% means that there is no threshold such that FDR is controlled at q.
%
%
% Background
%
% For a given threshold on a statistic image, the False Discovery Rate
% is the proportion of suprathreshold voxels which are false positives.
% Recall that the thresholding of each voxel consists of a hypothesis
% test, where the null hypothesis is rejected if the statistic is larger
% than threshold.  In this terminology, the FDR is the proportion of
% rejected tests where the null hypothesis is actually true.
%
% A FDR proceedure produces a threshold that controls the expected FDR
% at or below q.  The FDR adjusted p-value for a voxel is the smallest q
% such that the voxel would be suprathreshold.
%
% In comparison, a traditional multiple comparisons proceedure
% (e.g. Bonferroni or random field methods) controls Familywise Error
% rate (FWER) at or below alpha.  FWER is the *chance* of one or more
% false positives anywhere (whereas FDR is a *proportion* of false
% positives).  A FWER adjusted p-value for a voxel is the smallest alpha
% such that the voxel would be suprathreshold.
%
% If there is truely no signal in the image anywhere, then a FDR
% proceedure controls FWER, just as Bonferroni and random field methods
% do. (Precisely, controlling E(FDR) yeilds weak control of FWE).  If
% there *is* some signal in the image, a FDR method should be more powerful
% than a traditional method.
%
%
% References
%
% Benjamini & Hochberg (1995), "Controlling the False Discovery Rate: A
% Practical and Powerful Approach to Multiple Testing". J Royal Stat Soc,
% Ser B.  57:289-300.
%
% Benjamini & Yekutieli (2001), "The Control of the false discovery rate
% in multiple testing under dependency". To appear, Annals of Statistics.
% Available at http://www.math.tau.ac.il/~benja
%__________________________________________________________________________
% Copyright (C) 2002-2015 Wellcome Trust Centre for Neuroimaging

% Thomas Nichols
% $Id: spm_uc_FDR.m 6643 2015-12-11 16:57:25Z guillaume $


if (nargin<6), Vm = []; end

% Set Benjamini & Yeuketeli cV for independence/PosRegDep case
%--------------------------------------------------------------------------
cV = 1;


% Load, mask & sort statistic image (if needed)
%--------------------------------------------------------------------------
if isstruct(Vs)
    Ts = spm_data_read(Vs(1));
    for i=2:numel(Vs)
        Ts = min(Ts,spm_data_read(Vs(i)));
    end
    if ~isempty(Vm)
        if isstruct(Vm)
            Ts(spm_data_read(Vm)==0) = [];
        elseif numel(Vm)==1
            Ts(Ts==Vm) = [];
        else
            Ts = Ts(Vm);
        end
    end
    Ts(isnan(Ts)) = [];
    Ts = sort(Ts(:));
    if STAT ~= 'P', Ts = flipud(Ts); end
end


% Calculate p values of image (if needed)
%--------------------------------------------------------------------------
if isstruct(Vs)
    Ps = spm_z2p(Ts,df,STAT,n);
else
    Ps = Vs;
end

S = length(Ps);


% Calculate FDR inequality RHS
%--------------------------------------------------------------------------
Fi  = (1:S)'/S*q/cV;


% Find threshold
%--------------------------------------------------------------------------
I = find(Ps<=Fi, 1, 'last' );
if isempty(I)
    if STAT == 'P'
        u = 0;
    else
        u = Inf;
    end
else
    if isstruct(Vs)
        u = Ts(I);
    else
        % We don't have original statistic values; determine from p-value...
        if      STAT == 'Z'
            u = spm_invNcdf(1-Ps(I).^(1/n));
        elseif  STAT == 'T'
            u = spm_invTcdf(1-Ps(I).^(1/n),df(2));
        elseif  STAT == 'X'
            u = spm_invXcdf(1-Ps(I).^(1/n),df(2));
        elseif  STAT == 'F'
            u = spm_invFcdf(1-Ps(I).^(1/n),df);
            % ... except in case when we want a P-value
        elseif  STAT == 'P'
            u = Ps(I);
        end
    end
end
