function [dat,sett,mu] = spm_mb_fit(dat,sett)
% Multi-Brain - Groupwise normalisation and segmentation of images
% FORMAT [dat,sett,mu] = spm_mb_fit(dat,sett)
%
% OUTPUT
% dat                 - struct of length N storing each subject's information
% mu                  - array with template data
% sett  (inputParser) - struct storing final algorithm settings
% model (inputParser) - struct storing shape and appearance model
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


% Repeatable random numbers
%--------------------------------------------------------------------------
rng('default'); rng(1);

% If SPM has been compiled with OpenMP support then the number of threads
% are here set to speed up the algorithm
%--------------------------------------------------------------------------
if sett.nworker > 1
    setenv('SPM_NUM_THREADS', sprintf('%d',0));
else
    setenv('SPM_NUM_THREADS', sprintf('%d',-1));
end

% Delete any velocities and deformation files that could confuse things
%--------------------------------------------------------------------------
spm_mb_io('delete_files',dat)

% Get zoom (multi-scale) settings
%--------------------------------------------------------------------------
dmu     = sett.mu.d;
nz      = max(ceil(log2(min(dmu(dmu~=1))) - log2(sett.min_dim)), 1);
sz      = spm_mb_shape('zoom_settings', sett.v_settings, sett.mu, nz);
sett.ms = sz(end);

% Get template size and orientation
%--------------------------------------------------------------------------
te  = 0;
Mmu = sett.mu.Mmu;
if isfield(sett.mu,'exist')
    % Shrink given template
    mu0     = spm_mb_io('get_data', sett.mu.exist.mu);
    mu      = spm_mb_shape('shrink_template', mu0, Mmu, sett);
    updt_mu = false;
else
    % Random template
    nit_mu  = 1;
    mu      = bsxfun(@minus, randn([sett.ms.d sett.K], 'single'), ...
                             randn([sett.ms.d      1], 'single'));
    updt_mu = true;
end

nit_zm0   = 3;
updt_aff  = true;
updt_diff = all(isfinite(sett.v_settings));

% Specify subsampling
%--------------------------------------------------------------------------
for n=1:numel(dat)
    dat(n).samp  = get_samp(sett.ms.Mmu, dat(n).Mat, sett.sampdens);
    if isfield(dat(n).model, 'gmm')
        dat(n).model.gmm.samp = [1 1 1];
    end
end

% Update affine only
%--------------------------------------------------------------------------
fprintf('Rigid (zoom=1/%d): %d x %d x %d\n', 2^(numel(sz)-1), sett.ms.d);
spm_plot_convergence('Init', 'Rigid Alignment', 'Objective', 'Iteration');
E      = Inf;
for it0=1:128
    if it0>1
        oE  = E/nvox(dat);
    else
        oE  = Inf;
    end
    if updt_mu
        if it0<=12 && ~rem(it0,3), dat = spm_mb_appearance('restart', dat, sett); end
        [mu, sett, dat, te] = iterate_mean(mu, sett, dat, nit_mu);
    end

    % For visual debugging (disable/enable in debug_show_mu())
    debug_show_mu(mu, ['Affine (it=' num2str(it0) ' | N=' num2str(numel(dat)) ')']);

    % UPDATE: rigid
    dat   = spm_mb_shape('update_simple_affines',dat,mu,sett);
    sett  = spm_mb_appearance('update_prior',dat, sett);
    E     = sum(sum(cat(2,dat.E),2),1) + te;  % Cost function after previous update
    fprintf('%8.4f\n', E/nvox(dat));
    do_save(mu,sett,dat);
    spm_plot_convergence('Set',E/nvox(dat));

    % Finished rigid alignment?
    % Note that for limited field of view templates, the objective
    % function can increase as well as decrease.
    if it0>12 && abs(oE-E/nvox(dat)) < sett.tol
        countdown = countdown - 1;
        if countdown==0
            break;
        end
    else
        countdown = 6;
    end
end
spm_plot_convergence('Clear');


% Init shape model parameters
%--------------------------------------------------------------------------
dat    = spm_mb_shape('init_def',dat,sett);
nit_mu = 4;

% Update affine and diffeo (iteratively decreases the template resolution)
%--------------------------------------------------------------------------
spm_plot_convergence('Init','Diffeomorphic Alignment','Objective','Iteration');
for zm=numel(sz):-1:1 % loop over zoom levels
    fprintf('\nzoom=1/%d: %d x %d x %d\n', 2^(zm-1), sett.ms.d);

    % Specify subsampling
    for n=1:numel(dat)
        dat(n).samp  = [1 1 1];
        if isfield(dat(n).model, 'gmm')
            dat(n).model.gmm.samp = get_samp(sett.ms.Mmu, dat(n).Mat, sett.sampdens);
        end
    end

    if ~updt_mu
        mu = spm_mb_shape('shrink_template', mu0, Mmu, sett);
    else
        dat = spm_mb_appearance('restart', dat, sett);
        [mu, sett, dat, te] = iterate_mean(mu, sett, dat, nit_mu);
    end

    if updt_aff
        % UPDATE: rigid
        dat   = spm_mb_shape('update_affines', dat, mu, sett);
        sett  = spm_mb_appearance('update_prior', dat, sett);
        E     = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after previous update
        fprintf('%8.4f', E/nvox(dat));
        spm_plot_convergence('Set', E/nvox(dat));
    end
    fprintf('\n');

    nit_max = nit_zm0 + (zm - 1)*2;
    for it0=1:nit_max

        oE  = E/nvox(dat);
        if updt_mu
            [mu, sett, dat, te, E] = iterate_mean(mu, sett, dat, nit_mu);
        end

        % For visual debugging (disable/enable in debug_show_mu())
        debug_show_mu(mu, ['Diffeo (it=' num2str(zm) ', ' num2str(it0) ' | N=' num2str(numel(dat)) ')']);

        if updt_diff
            % UPDATE: diffeo
            dat   = spm_mb_shape('update_velocities', dat, mu, sett);
            sett  = spm_mb_appearance('update_prior', dat, sett);
            E     = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after previous update
            fprintf('%8.4f', E/nvox(dat));
            spm_plot_convergence('Set', E/nvox(dat));
        end
        fprintf('\n');

        if it0==nit_max || (oE-E/nvox(dat) < sett.tol && it0>=nit_zm0)
            break;
        else
            % Compute deformations from velocities (unless this is to be done
            % on the zoomed versions).
            if updt_diff
                dat   = spm_mb_shape('update_warps', dat, sett);
            end
        end
    end

    if zm > 1
        oMmu        = sett.ms.Mmu;
        d0          = size(mu);
        sett.ms     = copy_fields(sz(zm-1), sett.ms);
        if updt_diff
            dat     = spm_mb_shape('zoom_defs', dat, sett, oMmu, d0);
        end
        if updt_mu
            [mu,te] = spm_mb_shape('zoom_mean',  mu, sett, oMmu);
        end
    end

    if updt_diff
        dat = spm_mb_shape('update_warps', dat, sett); % Shoot new deformations
    end
    do_save(mu,sett,dat);
end
%spm_plot_convergence('Clear');
%==========================================================================

%==========================================================================
function nv = nvox(dat)
nv = sum([dat.nvox]);
%==========================================================================

%==========================================================================
function samp = get_samp(Mmu, Mf, sampdens)
if nargin<3
    n = 16;
else
    n = sampdens^3;
end
vmu   = sqrt(sum(Mmu(1:3,1:3).^2));
vf    = sqrt(sum( Mf(1:3,1:3).^2));
samp  = max(round(((prod(vmu)/prod(vf)/n).^(1/3))./vf),1);
samp  = min(samp,5);
%==========================================================================

%==========================================================================
function [mu, sett, dat, te, E] = iterate_mean(mu, sett, dat, nit_mu)
% UPDATE: mean
if nargin<4, nit_mu = 1; end
E      = Inf;
sampd  = mean_samp_dens(dat);
te     = spm_mb_shape('template_energy', mu, sett.ms.mu_settings, sampd);
for it=1:nit_mu
    [mu,dat] = spm_mb_shape('update_mean', dat, mu, sett, sampd);
    oE       = E;
    E        = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after previous update
    sett     = spm_mb_appearance('update_prior', dat, sett);
    te       = spm_mb_shape('template_energy', mu, sett.ms.mu_settings, sampd);
    fprintf('%8.4f', E/nvox(dat));
    spm_plot_convergence('Set',E/nvox(dat));
    do_save(mu,sett,dat);
    if it>1 && oE-E < sett.tol*nvox(dat); break; end
end
%==========================================================================

%==========================================================================
function s = mean_samp_dens(dat)
s = 0;
for n=1:numel(dat)
    s = s + prod(dat(n).samp);
end
s = s/numel(dat);
%==========================================================================

%==========================================================================
function to = copy_fields(from, to)
fn = fieldnames(from);
for i=1:numel(fn)
    to.(fn{i}) = from.(fn{i});
end
%==========================================================================

%==========================================================================
function do_save(mu, sett, dat)
if isfield(sett, 'save') && sett.save
    % Save results so far
    spm_mb_io('save_template', mu, sett);
    save(fullfile(sett.odir, ['mb_fit_' sett.onam '.mat']), 'sett', 'dat');
end
%==========================================================================

%==========================================================================
function debug_show_mu(mu, fig_title)
do_show = false;
spm_mb_appearance('debug_show', mu, 'template', 1, fig_title, do_show);
%==========================================================================
