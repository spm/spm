function spm_dcm_fmri_image(P)
% Image display of A, B, C and D coupling matrices
% FORMAT spm_dcm_fmri_image(P)
%
% P.A, P.B{1}, ...     - connections of weighted directed graph
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fmri_image.m 4627 2012-01-24 20:56:30Z karl $


% A and C - matrices
%--------------------------------------------------------------------------
subplot(6,2,1)
imagesc(P.A)
axis image
title('Average (A)','FontSize',16)
xlabel('from')
ylabel('to')

subplot(6,2,6 + 1)
bar(P.A)
title('Average (A)','FontSize',16)
xlabel('target')
spm_axis tight

subplot(6,2,2)
imagesc(P.C)
axis image
title('Exogenous (C)','FontSize',16)
xlabel('from')
ylabel('to')

subplot(6,2,6 + 2)
bar(P.C)
title('Exogenous (C)','FontSize',16)
xlabel('target')
spm_axis tight

% B - matrices
%--------------------------------------------------------------------------
m     = size(P.B,3);
for i = 1:m
    
    subplot(6,m,1*m + i)
    imagesc(P.B(:,:,i))
    axis image
    title('Bilinear (B)','FontSize',16)
    xlabel('from')
    ylabel('to')
    
    subplot(6,m,4*m + i)
    bar(P.B(:,:,i))
    title('Bilinear (B)','FontSize',16)
    xlabel('target')
    spm_axis tight
    
end

% D - matrices
%--------------------------------------------------------------------------
m     = size(P.D,3);
for i = 1:m
    subplot(6,m,2*m + i)
    imagesc(P.D(:,:,i))
    axis image
    title('Nonlinear (D)','FontSize',16)
    xlabel('from')
    ylabel('to')
    
    subplot(6,m,5*m + i)
    bar(P.D(:,:,i))
    title('Nonlinear (D)','FontSize',16)
    xlabel('target')
    spm_axis tight
end