function [Q] = spm_P_peakFDR(Z,df,STAT,R,n,ui,Ps)
% Return the corrected peak FDR q-value
% FORMAT [Q] = spm_P_peakFDR(Z,df,STAT,R,n,ui,Ps)
%
% Z        - height {minimum over n values}
% df       - [df{interest} df{residuals}]
% STAT     - Statistical field
%            'Z' - Gaussian field
%            'T' - T - field
%            'X' - Chi squared field
%            'F' - F - field
% R        - RESEL Count {defining search volume}
% n        - Conjunction number
% ui       - feature-inducing threshold
% Ps       - Vector of sorted (ascending) p-values
%
% Q        - FDR q-value
%__________________________________________________________________________
%
% References
% J.R. Chumbley and K.J. Friston, "False discovery rate revisited: FDR and 
% topological inference using Gaussian random fields". NeuroImage,
% 44(1):62-70, 2009.
%
% J.R. Chumbley, K.J. Worsley, G. Flandin and K.J. Friston, "Topological
% FDR for NeuroImaging". Under revision.
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Justin Chumbley & Guillaume Flandin
% $Id: spm_P_peakFDR.m 2764 2009-02-19 15:30:03Z guillaume $

% Expected Euler characteristic for level ui
%--------------------------------------------------------------------------
[P, p, Eu] = spm_P_RF(1, 0, ui, df, STAT, R, n);

% Expected Euler characteristic for level Z
%--------------------------------------------------------------------------
[P, p, Ez] = spm_P_RF(1, 0, Z, df, STAT, R, n);

% Uncorrected p-value for peaks using Random Field Theory
%--------------------------------------------------------------------------
Z = Ez / Eu;

% q value using the  Benjamini & Hochberch False Discovery Rate procedure
%--------------------------------------------------------------------------
Q = spm_P_FDR(Z, df, 'P',n, Ps);
