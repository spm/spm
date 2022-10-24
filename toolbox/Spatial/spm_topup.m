function VDM = spm_topup(vol1, vol2, FWHM, reg, outdir)
% Correct susceptibility distortions using topup
% FORMAT VDM = spm_topup(vol1, vol2, FWHM, reg, save)
% vol1   - path to first image (blip up/down image)
% vol2   - path to second image (blip down/up image)
% fwhm   - Gaussian kernel spatial scales (default: [8 4 2 1 0.1])
% reg    - regularisation settings (default: [0 10 100])
%          See spm_field for details:
%            - [1] Penalty on absolute values.
%            - [2] Penalty on the `membrane energy'. This penalises
%                  the sum of squares of the gradients of the values.
%            - [3] Penalty on the `bending energy'. This penalises
%                  the sum of squares of the 2nd derivatives.
% outdir - output directory
%
% VDM    - voxel displacement map
%
% Reference:
%
% J.L.R. Andersson, S. Skare, J. Ashburner. How to correct susceptibility
% distortions in spin-echo echo-planar images: application to diffusion
% tensor imaging. Neuroimage, 20(2):870-888, 2003.
% https://doi.org/10.1016/s1053-8119(03)00336-7
%__________________________________________________________________________

% John Ashburner, Nicole Labra Avila
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging


%-Optional input parameters
%--------------------------------------------------------------------------
if nargin < 3
    FWHM = [8 4 2 1 0.1]; % spatial scales
end
if nargin < 4
    reg = [0 10 100];     % regularisation
end
if nargin < 5
    outdir = '';          % output directory
end

%-Load images and estimate noise level
%--------------------------------------------------------------------------
P   = strvcat(fullfile(vol1), fullfile(vol2));
Nii = nifti(P);
vx  = sqrt(sum(Nii(1).mat(1:3,1:3).^2)); % Voxel sizes
reg = [vx reg];                          % Regularisation settings (spm_field)
sd  = spm_noise_estimate(P);             % Rice mixture model to estimate image noise
sig2 = sum(sd.^2);                       % Variance of difference is sum of variances

spm_field('bound',1);                    % Set boundary conditions for spm_field

% Set up identity transform (id) and displacements (u)
%--------------------------------------------------------------------------
d   = size(Nii(1).dat);
id  = zeros([d 3],'single');
[id(:,:,:,1),id(:,:,:,2),id(:,:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),single(1:d(3)));

u   = zeros(d,'single'); % Starting estimates of displacements
phi = id;                % x and z components are identity.
                         % y component will have displacement added/subtracted

%-Topup
%==========================================================================
for fwhm = FWHM % Loop over spatial scales
    fprintf('%3g:', fwhm);

    % Smooth the images according to the current spatial scale
    %----------------------------------------------------------------------
    f1  = single(Nii(1).dat(:,:,:));
    f2  = single(Nii(2).dat(:,:,:));
    spm_smooth(f1,f1,[1 1 1]*fwhm); % Note the side effects
    spm_smooth(f2,f2,[1 1 1]*fwhm); % Note the side effects

    n_acceptable = 0;
    E = 0;

    % Iterate until convergence (10 iterations max)
    %----------------------------------------------------------------------
    for it = 1:10

        % Sample the blip up/down image and its gradients in y direction
        phi(:,:,:,2) = id(:,:,:,2) + u;
        [wf1,~,d1,~] = spm_diffeo('bsplins',f1,phi,[1 1 1  0 0 0]);

        % Sample the blip down/up image and its gradients in y direction
        phi(:,:,:,2) = id(:,:,:,2) - u;
        [wf2,~,d2,~] = spm_diffeo('bsplins',f2,phi,[1 1 1  0 0 0]);

        % Regularisation term is \tfrac{1}{2} u^T L u. Compute L u. 
        gu  = spm_field('vel2mom', u, reg);

        % Compute cost function. Note the slightly ad hoc treatment of
        % missing data (where one of the phis points outside the FOV).
        % Note that there isn't a clear underlying generative model
        % of the data (so "log-likelihood" is in scare quotes).
        res = 0.5*((wf1 - wf2).^2/sig2 + gu.*u);
        msk = isfinite(res);
        Eprev = E;
        E   = sum(res(msk(:)))/sum(msk(:));
        fprintf(' %7.4f', E)

        % Gradient of "log-likelihood" term
        g   = (d1 + d2).*(wf1 - wf2)/sig2;
        
        % Diagonal of Hessian of "log-likelihood" term
        h   = (d1 + d2).^2/sig2;

        % Mask out missing data
        msk = ~isfinite(g);
        g(msk) = 0;
        h(msk) = 0;

        g   = g + gu; % Gradient of "log-likelihood" and log-prior terms

        % Gauss-Newton update
        % u \gets u - (diag(h) + L)\(g + L u)
        u   = u - spm_field(h, g, [reg 3 3]);

        % Check convergence
        if abs(E - Eprev) < 1e-3                % Stopping criterion
            n_acceptable = n_acceptable + 1;    % Closer to stopping
        else
            n_acceptable = 0;                   % Restart the countdown
        end
        if n_acceptable >= 3, break; end        % Done

        % Display topup iterations in Graphics window
        %------------------------------------------------------------------
        if true
            pl = ceil(size(f1,3)*0.2);
            FG = spm_figure('GetWin','Graphics');

            ax = subplot(2,2,1,'Parent',FG);
            imagesc(ax,[f1(:,:,pl)' f2(:,:,pl)']);
            axis(ax,'image','xy','off');
            title(ax,'Original');

            ax = subplot(2,2,2,'Parent',FG);
            imagesc(ax,[wf1(:,:,pl)' wf2(:,:,pl)']);
            axis(ax,'image','xy','off');
            title(ax,'Warped');

            ax = subplot(2,2,3,'Parent',FG);
            imagesc(ax,u(:,:,pl)');
            axis(ax,'image','xy','off');
            colorbar(ax);
            title(ax,'VDM');

            ax = subplot(2,2,4,'Parent',FG);
            imagesc(ax,g(:,:,pl)');
            axis(ax,'image','xy','off');
            colorbar(ax);
            title(ax,'dE/du');

            drawnow
        end
    end
    fprintf('\n');
end

%-Save distortion-corrected blip up/down and down/up images) and VDM
%==========================================================================

% wf1 (blip up/down image)
%--------------------------------------------------------------------------
basename = spm_file(vol1,'basename');
oname    = spm_file(basename,'prefix','w','ext','.nii');
oname    = fullfile(outdir,oname);
Nio      = nifti;                              
Nio.dat  = file_array(oname,d,'float32');    
Nio.mat  = Nii(1).mat;          
Nio.mat0 = Nii(1).mat;
create(Nio);         
Nio.dat(:,:,:) = wf1; 

% wf2 (blip down/up image)
%--------------------------------------------------------------------------
basename = spm_file(vol2,'basename');
oname    = spm_file(basename,'prefix','w','ext','.nii');
oname    = fullfile(outdir,oname);
Nio      = nifti; 
Nio.dat  = file_array(oname,d,'float32'); 
Nio.mat  = Nii(2).mat; 
Nio.mat0 = Nii(2).mat0;
create(Nio);          
Nio.dat(:,:,:) = wf2;

% Write VDM file
%--------------------------------------------------------------------------
basename       = spm_file(vol1,'basename');
oname          = spm_file(basename,'prefix','vdm5_','ext','.nii');
oname          = fullfile(outdir,oname);
Nio            = nifti;  
Nio.dat        = file_array(oname,size(u),'float32');
Nio.mat        = Nii(1).mat;  
create(Nio);          
Nio.dat(:,:,:) = u; 

VDM = Nio;
