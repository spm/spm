function [g] = spm_gx_fnirs(x,u,P,M)
% fNIRS optics equation
% FORMAT [g] = spm_gx_fnirs(x,u,P,M)
%
% x     state vector     (see spm_fx_fnirs)
% u     experimental inputs
% P     prior of latent variables
% M    model structure
%
% g     optical density changes
%___________________________________________________________________________
% References for optics equations:
% 1. Arridge, SR 1999. Optical tomography in medical imaging. Inverse Prob.
% 15: R41-R93.
% 2. Gagnon L, Yucel, MA, Dehaes, M, Cooper, RJ, Perdue, KL, Selb, J, Huppert TJ,
% Hoge RD, Boas DA, 2012. Quantification of the cortical contribution to
% the NIRS signal over NIRS-fMRI measurements. NeuroImage 59: 3933-3940.
% 3. Tak, S, Kempny, AM, Friston, KJ, Leff, AP, Penny WD, Dynamic causal
% modelling for functional near-infrared spectroscopy. NeuroImage 111: 338-349.
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny & Sungho Tak
% $Id: spm_gx_fnirs.m 6754 2016-03-25 06:44:58Z will $

% exponentiation of hemodynamic state variables
x(:, 3:6) = exp(x(:,3:6));

% optics parameters
%------------------------------------------------------------------------------------
%   N(1) - oxygen saturation                                            SO2
%   N(2) - baseline total-hb concentration [mM]                P0
%   N(3) - cortical weighting factor                                   w
%------------------------------------------------------------------------------------
N = [0.65 0.071 2];

% baseline dxy-hb concentration
base_hbr = N(2) .* (1-N(1));

% calculate changes in hemoglobin concentration
dq = (x(:,5)-1).*base_hbr; % dxy-Hb
dp = (x(:,6)-1).*N(2); % total Hb
dh = dp - dq; % oxy-Hb

% parameter for correction of pial vein contamination
if size(P.pv, 1) == (M.nch*2)
    % additionally correct pial vein oxygenation
    pv = N(3) .* exp(P.pv);
    pv = reshape(pv, M.nch, 2);
else
    pv = N(3).* exp(P.pv);
    pv = kron(ones(1,2), pv);
end

% calculate absorption coefficients
if M.rs > 0 % spatially distibuted source
    sh = []; sq = [];
    ns = size(x, 1); % number of sources
    for i = 1:ns
        rs = M.rs .* exp(P.rs(i));
        indx = find(M.A(i).dist_h < rs);
        
        s = (rs/ sqrt(8 * log(2))).^2 + eps;
        skrn = zeros(length(M.A(i).dist_h),1); % smoothing kernel
        skrn(indx) = exp(-(M.A(i).dist_h(indx).^2)/(2*s));
        
        sh = [sh; skrn * dh(i,1)];
        sq = [sq; skrn * dq(i,1)];
    end
else % point source
    sh = dh; sq = dq;
end

% calculate optical density changes
for i = 1:M.nwav
    S = spm_vec(M.A.(sprintf('S%d', i)));
    S = reshape(S, M.nch, []);
    g(:,i) = (pv .* (S * [sh sq])) * M.acoef(i,:)';
end
g = spm_vec(g);

