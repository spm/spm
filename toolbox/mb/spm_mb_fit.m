function [dat,sett,mu] = spm_mb_fit(dat,sett)
% Multi-Brain - Groupwise normalisation and segmentation of images
%
% FORMAT [dat,mu,sett] = spm_mb_fit(dat,sett)
%
% OUTPUT
% dat                 - struct of length N storing each subject's information
% mu                  - array with template data
% sett  (inputParser) - struct storing final algorithm settings
% model (inputParser) - struct storing shape and appearance model
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% $Id: spm_mb_fit.m 7907 2020-07-23 16:10:52Z john $


% Repeatable random numbers
%--------------------------------------------------------------------------
rng('default'); rng(1);

% If SPM has been compiled with OpenMP support then the number of threads
% are here set to speed up the algorithm
%--------------------------------------------------------------------------
if sett.nworker > 1
    setenv('SPM_NUM_THREADS',sprintf('%d',0));
else
    setenv('SPM_NUM_THREADS',sprintf('%d',-1));
end

% Get template size and orientation
%--------------------------------------------------------------------------
if isfield(sett.mu,'exist')
    mu0 = spm_mb_io('get_data',sett.mu.exist.mu);
end
Mmu = sett.mu.Mmu;

% Get zoom (multi-scale) settings
%--------------------------------------------------------------------------
dmu     = sett.mu.d;
nz      = max(ceil(log2(min(dmu(dmu~=1))) - log2(sett.min_dim)),1);
sz      = spm_mb_shape('zoom_settings',sett.v_settings,sett.mu,nz);
sett.ms = sz(end);

% Init shape model parameters
%--------------------------------------------------------------------------
dat = spm_mb_shape('init_def',dat,sett.ms);

% Init template
%--------------------------------------------------------------------------
if exist('mu0','var')
    % Shrink given template
    mu = spm_mb_shape('shrink_template',mu0,Mmu,sett);
else
    % Random template
    mu = randn([sett.ms.d sett.K],'single')*1.0;
end

% Start algorithm
%--------------------------------------------------------------------------
nit_zm0   = 3;
nit_aff   = 128;
updt_aff  = true;
updt_diff = all(isfinite(sett.v_settings));
updt_mu   = ~exist('mu0','var');
if ~exist('mu0','var')
    te = spm_mb_shape('template_energy',mu,sett.ms.mu_settings);
    E  = [Inf Inf];
else
    te = 0;
    E  = Inf;
end

% Update affine only
% No update of intensity priors during this phase.
%--------------------------------------------------------------------------
for n=1:numel(dat)
    dat(n).samp  = get_samp(sett.ms.Mmu,dat(n).Mat,sett.sampdens);
    if isfield(dat(n).model,'gmm')
        dat(n).model.gmm.samp = [1 1 1];
    end
end
updt_int = 'update_prior';
fprintf('Rigid (zoom=%d): %d x %d x %d\n',2^(numel(sz)-1),sett.ms.d);
spm_plot_convergence('Init','Rigid Alignment','Objective','Iteration');
E      = Inf;
for it0=1:nit_aff
    if it0>1
        oE  = E/nvox(dat);
    else
        oE  = Inf;
    end
    if updt_mu
        [mu,sett,dat,te,E] = iterate_mean(mu,sett,dat,te,E,updt_int);
    end

    if true
        % UPDATE: rigid
        dat   = spm_mb_shape('update_simple_affines',dat,mu,sett);
        E     = sum(sum(cat(2,dat.E),2),1) + te;  % Cost function after previous update
        sett  = spm_mb_appearance(updt_int,dat, sett);

        fprintf('%12.4e', E/nvox(dat));
    end
    fprintf('\n');
    do_save(mu,sett,dat);
    spm_plot_convergence('Set',E/nvox(dat));

    % Finished rigid alignment?
    % Note that for limited field of view templates, the objective
    % function can increase as well as decrease.
    if abs(oE-E/nvox(dat)) < sett.tol
        countdown = countdown - 1;
        if countdown==0
            break;
        end
    else
        countdown = 6;
    end
end
spm_plot_convergence('Clear');


% Update affine and diffeo (iteratively decreases the template resolution)
%--------------------------------------------------------------------------
spm_plot_convergence('Init','Diffeomorphic Alignment','Objective','Iteration');
for zm=numel(sz):-1:1 % loop over zoom levels
    fprintf('\nzoom=%d: %d x %d x %d\n', 2^(zm-1), sett.ms.d);
   %spm_plot_convergence('Init',['Diffeomorphic Alignment (' num2str(2^(zm-1)) ')'],'Objective','Iteration');
    for n=1:numel(dat)
        dat(n).samp  = [1 1 1];
        if isfield(dat(n).model,'gmm')
            dat(n).model.gmm.samp = get_samp(sett.ms.Mmu,dat(n).Mat,sett.sampdens);
        end
    end

    if ~updt_mu
        mu = spm_mb_shape('shrink_template',mu0,Mmu,sett);
    else
        [mu,sett,dat,te,E] = iterate_mean(mu,sett,dat,te,E);
    end

    if updt_aff
        % UPDATE: rigid
        dat   = spm_mb_shape('update_affines',dat,mu,sett);
        E     = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after previous update
        sett  = spm_mb_appearance('update_prior',dat, sett);
        fprintf('%12.4e', E/nvox(dat));
        spm_plot_convergence('Set',E/nvox(dat));
    end
    fprintf('\n');

    EE     = inf(1,sum([updt_mu&&updt_diff, updt_diff])); % For tracking objfun
    nit_zm = nit_zm0 + (zm - 1); % use nit_zm0 only for zm = 1
    for it0=1:nit_zm

       %oEE = EE;
        i   = 1;   % For tracking objfun

        if updt_mu
            [mu,sett,dat,te,E] = iterate_mean(mu,sett,dat,te,E);
            EE(i)  = E;
            i      = i+1;
        end

        if updt_diff
            % UPDATE: diffeo
            dat   = spm_mb_shape('update_velocities',dat,mu,sett);
            E     = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after previous update
            sett  = spm_mb_appearance('update_prior',dat, sett);
            EE(i) = E;
            fprintf('%12.4e', E/nvox(dat));
            spm_plot_convergence('Set',E/nvox(dat));
        end
        fprintf('\n');

        % Compute deformations from velocities (unless this is to be done
        % on the zoomed versions).
        if it0<nit_zm && updt_diff
            dat   = spm_mb_shape('update_warps',dat,sett);
        end
    end

    if zm > 1
        oMmu           = sett.ms.Mmu;
        sett.ms        = copy_fields(sz(zm-1), sett.ms);
        if      updt_mu && updt_diff
            [dat,mu]   = spm_mb_shape('zoom_volumes',dat,mu,sett,oMmu);
        elseif  updt_mu && ~updt_diff
            [~,mu]     = spm_mb_shape('zoom_volumes',[],mu,sett,oMmu);
        elseif ~updt_mu && updt_diff
            dat        = spm_mb_shape('zoom_volumes',dat,mu,sett,oMmu);
        end
        if updt_mu, te = spm_mb_shape('template_energy',mu,sett.ms.mu_settings); end % Compute template energy
    end

    if updt_diff
        dat        = spm_mb_shape('update_warps',dat,sett); % Shoot new deformations
    end
    do_save(mu,sett,dat);
end
%spm_plot_convergence('Clear');
%==========================================================================

%==========================================================================
function nv = nvox(dat)
nv = 0;
for n=1:numel(dat)
    nv = nv+dat(n).nvox;
end
%==========================================================================

%==========================================================================
function samp = get_samp(Mmu,Mf,sampdens)
if nargin<3
    n=16;
else
    n = sampdens^3;
end
vmu  = sqrt(sum(Mmu(1:3,1:3).^2));
vf   = sqrt(sum( Mf(1:3,1:3).^2));
samp = max(round(((prod(vmu)/prod(vf)/n).^(1/3))./vf),1);
samp = min(samp,5);
%==========================================================================

%==========================================================================
function [mu,sett,dat,te,E] = iterate_mean(mu,sett,dat,te,E,updt_int)
% UPDATE: mean
if nargin<6, updt_int = 'update_prior'; end
nit_mu = 5;
for it=1:nit_mu
    [mu,dat] = spm_mb_shape('update_mean',dat, mu, sett);
    oE       = E;
    E        = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after previous update
    sett     = spm_mb_appearance(updt_int,dat, sett);
    te       = spm_mb_shape('template_energy',mu,sett.ms.mu_settings);
    fprintf('%12.4e', E/nvox(dat));
    spm_plot_convergence('Set',E/nvox(dat));
    do_save(mu,sett,dat);
    if it>1 && oE-E < sett.tol*nvox(dat); break; end
end
%==========================================================================

%==========================================================================
function to = copy_fields(from,to)
fn = fieldnames(from);
for i=1:numel(fn)
    to.(fn{i}) = from.(fn{i});
end
%==========================================================================

%==========================================================================
function do_save(mu,sett,dat)
if isfield(sett,'save') && sett.save
    % Save results so far
    spm_mb_io('save_template',mu,sett);
   %sett = rmfield(sett,{'ms'});
    save(fullfile(sett.odir,['mb_fit_' sett.onam '.mat']),'sett','dat');
end
%==========================================================================
