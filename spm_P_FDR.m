function [P,Ps,Ts] = spm_P_FDR(Z,df,STAT,n,Vs,Vm)
% Returns the corrected FDR P value 
% FORMAT [P,Ps,Ts] = spm_P_FDR(Z,df,STAT,n,Vs[,Vm])
%
% Z     - height (minimum of n statistics)
% df    - [df{interest} df{error}]
% STAT  - Statisical field
%		'Z' - Gaussian feild
%		'T' - T - feild
%		'X' - Chi squared feild
%		'F' - F - feild
% n     - number of component SPMs in conjunction
% Vs    - Mapped statistic image(s)
%          -or-
%         Vector of sorted p-values (saves i/o w/ repeated calls)
% Vm    - Mask in 1 of 3 forms
%           o Scalar, indicating implicit mask value in statistic image(s)
%           o Vector of indicies of elements within mask
%           o Mapped mask image
%
%
% P     - corrected FDR   P value
% Ps    - Sorted p-values
% Ts    - Sorted statistic values
%
%___________________________________________________________________________
%
% The Benjamini & Hochberch (1995) False Discovery Rate (FDR) procedure
% finds a threshold u such that the expected FDR is at most q.  spm_P_FDR
% returns the smallest q such that Z>u. 
%
% For repeated use of a given statistic image, return Ps in the place
% of Vs:
%      [P Ps] = spm_P_FDR(Z1,df,STAT,n,Vs);  %-Initialization, image read
%      P      = spm_P_FDR(Z2,df,STAT,n,Ps);  %-Image not read, Ps used
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
% rate (FWER) at or below alpha.  FWER is the chance of one or more
% false positives *anywhere* (not just among suprathreshold voxels).  A
% FWER adjusted p-value for a voxel is the smallest alpha such that the
% voxel would be suprathreshold.
%
% 
% If there is truely no signal in the image anywhere, then a FDR
% proceedure controls FWER, just as Bonferroni and random field methods
% do. (Precisely, controlling E(FDR) yeilds weak control of FWE).  If
% there *is* some signal in the image, a FDR method will be more powerful
% than a traditional method.
%
%
% For careful definition of FDR-adjusted p-values (and distinction between
% corrected and adjusted p-values) see Yekutieli & Benjamini (1999).
%
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
%
% Yekutieli & Benjamini (1999). "Resampling-based false discovery rate
% controlling multiple test proceedures for correleated test
% statistics".  J of Statistical Planning and Inference, 82:171-196.
%___________________________________________________________________________
% %W% Thomas Nichols %E%

if (nargin<6), Vm = []; end

% Set Benjamini & Yeuketeli cV for independence/PosRegDep case
%-----------------------------------------------------------------------
cV = 1; 


% Load, mask & sort statistic image (if needed)
%-----------------------------------------------------------------------
if isstruct(Vs)
  if (n ~= length(Vs))
    error(sprintf('n & number of mapped images doesn''t match (%d,%d)',...
		  n,length(Vs)));
  end
  Ts = spm_read_vols(Vs(1));
  for i=2:n
    Ts = min(Ts,spm_read_vols(Vs(i)));
  end
  Ts(isnan(Ts)) = [];
  if ~isempty(Vm)
    if isstruct(Vm)
      Ts(spm_read_vols(Vm)==0) = [];
    elseif (prod(size(Vm))==1)
      Ts(Ts==Vm) = [];
    else 
      Ts = Ts(Vm);
    end
  end
  Ts = flipud(sort(Ts(:)));
end


% Calculate p values of image (if needed)
%-----------------------------------------------------------------------
if isstruct(Vs)
  if      STAT == 'Z'
    Ps = (1-spm_Ncdf(Ts)).^n;
  elseif  STAT == 'T'
    Ps = (1-spm_Tcdf(Ts,df(2))).^n;
  elseif  STAT == 'X'
    Ps = (1-spm_Xcdf(Ts,df(2))).^n;
  elseif  STAT == 'F'
    Ps = (1-spm_Fcdf(Ts,df)).^n;
  end
else
  Ps = Vs;
end

S = length(Ps);


% Calculate p value of Z
%-----------------------------------------------------------------------
if      STAT == 'Z'
  PZ = (1-spm_Ncdf(Z)).^n;
elseif  STAT == 'T'
  PZ = (1-spm_Tcdf(Z,df(2))).^n;
elseif  STAT == 'X'
  PZ = (1-spm_Xcdf(Z,df(2))).^n;
elseif  STAT == 'F'
  PZ = (1-spm_Fcdf(Z,df)).^n;
end


%-Calculate FDR p values
%
% If Z is a value in the statistic image, then the adjusted p-value
% defined in Yekutieli & Benjamini (1999) (eqn 3) is obtained.  If Z
% isn't a value in the image, then the adjusted p-value for the next
% smallest statistic value (next largest uncorrected p) is returned.
%
%-----------------------------------------------------------------------

Qs  = Ps*S./(1:S)'*cV;    %-"Corrected" p-values

P   = zeros(size(Z));
for i=1:length(Z)

  % Find PZ(i) in Ps, or smallest Ps(j) s.t. Ps(j) >= PZ(i)
  I = min(find(Ps>=PZ(i)));  

  % Find "adjusted" p-values
  if isempty(I)
    P(i) = 1;
  else
    P(i) = min(Qs(I:S));  %-"Adjusted" p-values
  end

end
