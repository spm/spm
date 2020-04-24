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
%
%__________________________________________________________________________
%
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging

% $Id$

% Repeatable random numbers
rng('default'); rng(1);

%------------------
% Get template size and orientation
%------------------
if isfield(sett.mu,'exist')
    mu0 = spm_mb_io('GetData',sett.mu.exist.mu);
end
Mmu = sett.mu.Mmu;

%------------------
% Get zoom (multi-scale) settings
%------------------
dmu     = sett.mu.d;
nz      = max(ceil(log2(min(dmu(dmu~=1))) - log2(sett.min_dim)),1);
sz      = spm_mb_shape('ZoomSettings',sett.v_settings,sett.mu,nz);
sett.ms = sz(end);

%------------------
% Init shape model parameters
%------------------
dat = spm_mb_shape('InitDef',dat,sett.ms);

%------------------
% Init template
%------------------
if exist('mu0','var')
    % Shrink given template
    mu = spm_mb_shape('ShrinkTemplate',mu0,Mmu,sett);
else
    % Random template
    mu = randn([sett.ms.d sett.K],'single')*1.0;
end

%------------------
% Start algorithm
%------------------
nit_zm0    = 3;
nit_aff    = 128;
updt_aff   = true;
updt_diff  = all(isfinite(sett.v_settings));
updt_mu    = ~exist('mu0','var');
if ~exist('mu0','var')
    te = spm_mb_shape('TemplateEnergy',mu,sett.ms.mu_settings);
    E  = [Inf Inf];
else
    te = 0;
    E  = Inf;
end

%------------------
% Update affine only
% No update of intensity priors during this phase.
%------------------
sett.tol = sett.tol*0.1;
for n=1:numel(dat)
    dat(n).samp  = GetSamp(sett.ms.Mmu,dat(n).Mat,sett.sampdens);
    dat(n).samp2 = [1 1 1];
end
updt_int = 'UpdatePrior';
fprintf('Rigid (zoom=%d): %d x %d x %d\n',2^(numel(sz)-1),sett.ms.d);
spm_plot_convergence('Init','Rigid Alignment & Burn In','Objective','Iteration');
EE     = inf(1,sum([updt_mu, 1])); % For tracking objfun
for it0=1:nit_aff
    oEE = EE;
    i   = 1;   % For tracking objfun

    if updt_mu
        [mu,sett,dat,te,E] = IterateMean(mu,sett,dat,te,E,updt_int);
        EE(i) = E;
        i     = i + 1;
    end

    if true
        % UPDATE: rigid
        dat   = spm_mb_shape('UpdateSimpleAffines',dat,mu,sett);
        E     = sum(sum(cat(2,dat.E),2),1) + te;  % Cost function after previous update
        sett  = spm_mb_appearance(updt_int,dat, sett);

        EE(i) = E;
        fprintf('%13.5e', E);
    end
    fprintf('\n');
    do_save(mu,sett,dat);
    spm_plot_convergence('Set',E);

    % Check convergence
    change = mean(abs(oEE - EE)./abs(EE));

    % Finished rigid alignment?
    if change < sett.tol
        countdown = countdown - 1;
        if countdown==0
            break;
        end
    else
        countdown = 4;
    end
end

sett.tol = sett.tol*10;
spm_plot_convergence('Clear');


%------------------
% Update affine and diffeo (iteratively decreases the template resolution)
%------------------
for zm=numel(sz):-1:1 % loop over zoom levels
    fprintf('\nzoom=%d: %d x %d x %d\n', 2^(zm-1), sett.ms.d);
    spm_plot_convergence('Init',['Diffeomorphic Alignment (' num2str(2^(zm-1)) ')'],'Objective','Iteration');
    for n=1:numel(dat),
        dat(n).samp  = [1 1 1];
        dat(n).samp2 = GetSamp(sett.ms.Mmu,dat(n).Mat,sett.sampdens);
    end

    if ~updt_mu
        mu = spm_mb_shape('ShrinkTemplate',mu0,Mmu,sett);
    else
        [mu,sett,dat,te,E] = IterateMean(mu,sett,dat,te,E);
    end

    if updt_aff
        % UPDATE: rigid
        dat   = spm_mb_shape('UpdateAffines',dat,mu,sett);
        E     = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after previous update
        sett  = spm_mb_appearance('UpdatePrior',dat, sett);
        fprintf('%13.5e', E);
        spm_plot_convergence('Set',E);
    end
    fprintf('\n');

    EE     = inf(1,sum([updt_mu&&updt_diff, updt_diff])); % For tracking objfun
    nit_zm = nit_zm0 + (zm - 1); % use nit_zm0 only for zm = 1
    for it0=1:nit_zm

        oEE = EE;
        i   = 1;   % For tracking objfun

        if updt_mu
            [mu,sett,dat,te,E] = IterateMean(mu,sett,dat,te,E);
            EE(i)  = E;
            i      = i+1;
        end

        if updt_diff
            % UPDATE: diffeo
            dat   = spm_mb_shape('UpdateVelocities',dat,mu,sett);
            E     = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after previous update
            sett  = spm_mb_appearance('UpdatePrior',dat, sett);
            EE(i) = E;
            fprintf('%13.5e', E);
            spm_plot_convergence('Set',E);
        end
        fprintf('\n');

       %% Check convergence and terminate if done
       %change = mean(abs(oEE - EE)./abs(EE));        
       %if change < sett.tol, break; end

        % Compute deformations from velocities (unless this is to be done
        % on the zoomed versions).
        if it0<nit_zm && updt_diff
            dat   = spm_mb_shape('UpdateWarps',dat,sett);
        end
    end

    if zm > 1
        oMmu           = sett.ms.Mmu;
        sett.ms        = CopyFields(sz(zm-1), sett.ms);
        if      updt_mu && updt_diff
            [dat,mu]   = spm_mb_shape('ZoomVolumes',dat,mu,sett,oMmu);
        elseif  updt_mu && ~updt_diff
            [~,mu]     = spm_mb_shape('ZoomVolumes',[],mu,sett,oMmu);
        elseif ~updt_mu && updt_diff
            dat        = spm_mb_shape('ZoomVolumes',dat,mu,sett,oMmu);
        end
        if updt_mu, te = spm_mb_shape('TemplateEnergy',mu,sett.ms.mu_settings); end % Compute template energy
    end

    if updt_diff
        dat        = spm_mb_shape('UpdateWarps',dat,sett); % Shoot new deformations
    end
    do_save(mu,sett,dat);
    spm_plot_convergence('Clear');
end
end
%==========================================================================

%==========================================================================
function samp = GetSamp(Mmu,Mf,sampdens)
if nargin<3, n=16; else n = sampdens^3; end
vmu  = sqrt(sum(Mmu(1:3,1:3).^2));
vf   = sqrt(sum( Mf(1:3,1:3).^2));
samp = max(round(((prod(vmu)/prod(vf)/n).^(1/3))./vf),1);
samp = min(samp,5);
end
%==========================================================================

%==========================================================================
function [mu,sett,dat,te,E] = IterateMean(mu,sett,dat,te,E,updt_int)
% UPDATE: mean
if nargin<6, updt_int = 'UpdatePrior'; end
nit_mu = 5;
for it=1:nit_mu
    [mu,dat] = spm_mb_shape('UpdateMean',dat, mu, sett);
    oE       = E;
    E        = sum(sum(cat(2,dat.E),2),1) + te; % Cost function after previous update
    sett     = spm_mb_appearance(updt_int,dat, sett);
    te       = spm_mb_shape('TemplateEnergy',mu,sett.ms.mu_settings);
    fprintf('%13.5e', E);
    spm_plot_convergence('Set',E);
    if it>1 && (oE-E)/abs(E) < sett.tol; break; end
end
end
%==========================================================================

%==========================================================================
function to = CopyFields(from,to)
fn = fieldnames(from);
for i=1:numel(fn)
    to.(fn{i}) = from.(fn{i});
end
end
%==========================================================================

%==========================================================================
function do_save(mu,sett,dat)
if isfield(sett,'save') && sett.save
    % Save results so far
    spm_mb_io('SaveTemplate',mu,sett);
   %dat  = rmfield(dat,{'samp','samp2'});
   %sett = rmfield(sett,{'ms'});
    save(fullfile(sett.odir,['mb_fit_' sett.onam '.mat']),'sett','dat');
end
end
%==========================================================================
