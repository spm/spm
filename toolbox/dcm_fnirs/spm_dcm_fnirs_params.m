function [A, B, C] = spm_dcm_fnirs_params(DCM)
% Calculate DCM parameters using estimated latent variables
% FORMAT [A, B, C] = spm_dcm_fnirs_params(DCM)
%
% DCM  - DCM structure (see spm_dcm_ui)
%
% A - Endogenous (fixed) connections
% B - Connections modulated by input
% C - Influence of input on regional activity 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny & Sungho Tak
% $Id: spm_dcm_fnirs_params.m 6754 2016-03-25 06:44:58Z will $

% intrinsic connection in the absence of input
r = 10.^5;

mlA = diag(DCM.Ep.A);
slA = sqrt(diag(DCM.Vp.A));

n = DCM.n;
for i = 1:n,
    latentA(i,:) = mlA(i) + slA(i) .* randn(1,r);
end

iA = (-0.5) .* exp(latentA); % intrinsic connectivity
miA = mean(iA, 2);
siA = std(iA, 0, 2);

nc = size(DCM.Ep.B, 3); % experimental conditions
miB = zeros(n, nc);
siB = zeros(n, nc);

for i = 1:nc
    mlB = diag(DCM.Ep.B(:,:, i)); % latent variables
    slB = sqrt(diag(DCM.Vp.B(:,:,i)));
    for j =1:n
        latentB = mlB(j) + slB(j) .* randn(1,r);
        iB = (-0.5) .* exp(latentA(j,:) + latentB) + 0.5 .* exp(latentA(j,:));
        miB(j,i) = mean(iB);
        siB(j,i) = std(iB);
    end
end

mA = DCM.Ep.A;
sA = sqrt(DCM.Vp.A);

indx = find(eye(n));
mA(indx) = miA;
sA(indx) = siA;

mB = reshape(DCM.Ep.B, n^2, nc)';
sB = reshape(sqrt(DCM.Vp.B), n^2, nc)';
for i = 1:nc,
    mB(i, indx) = miB(:,i)';
    sB(i, indx) = siB(:,i)';
end
mB = reshape(mB', n, n, nc);
sB = reshape(sB', n, n, nc);

mC = DCM.Ep.C./16;
sC = sqrt(DCM.Vp.C)./16;

A.mean = mA; 
A.std = sA; 

B.mean = mB; 
B.std = sB; 

C.mean = mC; 
C.std = sC; 



