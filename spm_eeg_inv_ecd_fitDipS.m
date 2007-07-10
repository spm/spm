function [sdip,fit_opt] = spm_eeg_inv_ecd_fitDipS(V,model,Vbr,n_dip,n_seeds,fit_opt,disp_f,MCS_f)

% FORMAT [sdip,set_loc_o] = spm_eeg_inv_ecd_fitDipS(V,model,Vbr,n_dip,n_seeds,fit_opt,disp_f,MCS_f)
%
% Fit  moving dipoles on EEG data
% Input :
%   - V:     EEG data (n_el x n_tbin)
%   - model: head model structure with: tessalated surfaces, electrodes setup,
%					conductivity values, spheres defintion.
%   - Vbr:   brain mask to limit the dipoles location.
%   - n_dip: nr of dipoles to be fitted.
%   - n_seeds: nr of random starting locations to be considered.
%   - fit_opt: option for the optimization
%       + q_model: model used: realistic BEM (1) or analytical 3 fitted sphere (2)
%       + or_opt: ways of fixing (or not) dipoles orientation over time
%                    (1: free orientation, 2: fixed or.as weighted mean orientation,
%                        3: fixed or as at EEG max power, 4: fixed a priori)
%       + rem_electr: electrodes to remove from the model (electrodes structure and inthe IFS matrices)
%       + q_fxd_loc: fix the dipoles onto their a priori location (1) or not (0)
%       + set_loc_o: priors for the dipole location
%       + q_fxd_or: fix the dipoles orientation as specified (1) or not (0)
%       + fxd_or: fixed orientation of the dipoles (used with or_opt=4)
%   - disp_f: display message after fitting (1) or not (0) (default 1)
%   - MCS_f: flag for Monte Carlo simulations, MCS, (1) or not (0). (default 0)
%
% Output :
%   - sdip: fitted dipole(s) strucure
%       + n_seeds: # seeds set used
%       + n_dip: # fitted dipoles on the EEG time series
%       + loc: location of fitted dipoles, cell{1,n_seeds}(3 x n_dip)
%       + L: leadfiled associated with each set of n_dip fitted dipoles, 
%            cell{1,n_seeds}(N_el x 3*n_dip)
%       + j: sources amplitude over the time window, cell{1,n_seeds}(3*n_dip x Ntimebins)
%       + res: residuals of each fit, vect(1,n_seeds)
%       + cost: cost of each fit, vect(1,n_seeds)
%       + rres: relative residuals of each fit, vect(1,n_seeds)
%       + varexpl: variance explained by fitted dipole(s)
%       + M: 4x4 affine transformation matrix mapping from brain image voxel coordinates
%            to real world coordinates
%       + Mtb: index of maximum power in EEG time series used
%       + fit_opt: fitting options used
%       + wind_be: time window defining timeseries used, in complete data set
%       + Pdata: data file
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips,
% $Id$


global MODEL V_BR V_EEG OR_OPT FXD_OR sdip
% Use sdip as global for the Monte Carlo tests.

if nargin<8
    MCS_f = 0;
end

if nargin<7
    disp_f = 1;
end

if nargin<6,
    fit_opt.q_model = 2;
    fit_opt.or_opt = 1;
    fit_opt.rem_electr = [];
    fit_opt.q_fxd_loc = 0;
    fit_opt.set_loc_o = [];
    fit_opt.q_fxd_or = 0;
    fit_opt.fxd_or = [];
end
if nargin<5, n_seeds = 1; end
if nargin<4, n_dip = 1; end
if nargin<3
    error('Provide at least, 1. EEG data, 2. the head model & 3. brain mask')
end
if ~isstruct(Vbr), Vbr = spm_vol(Vbr); end


% Ensure repeatable solutions as function of localisation number (ie seed) (RH)
%------------------------------------------------------------------------------
try
    val = D.val;
catch
    val = 1;
    D.val = val;
end
rand('state',val)

% Take average reference for EEG data (RH)
%-----------------------------------------
% (since leadfield mean-corrected in spm_eeg_inv_ecd_costSd.m)

V = V - repmat(mean(V,1),size(V,1),1);

% Pass the global vars
%---------------------
V_BR = Vbr;
V_EEG = V;
OR_OPT = fit_opt.or_opt;
FXD_OR = fit_opt.fxd_or;

% Generate random seeds, if a priori location are not provided...
%-------------------------------------------------------------
if ~MCS_f
    if isempty(fit_opt.set_loc_o)
        % if only 1 seed & 1 dipole, let it at [0 0 0]'mm.
        % if only 1 seed & x dipoles, take them at symetric locations I specify.
        set_loc_1seed = [[-30 0 0]' [30 0 0]' [0 0 30]' [0 30 0]' [0 -30 0]'  [0 0 0]'];
        % This will be ok for up to 6 dipoles
        % Otherwise pick at random in brain volume.
        if n_seeds==1 && n_dip==1
            fit_opt.set_loc_o = zeros(3,1);
        elseif n_seeds==1 && n_dip>1
            fit_opt.set_loc_o = set_loc_1seed(:,1:n_dip);
        else
            loc_rand = rand(n_seeds*20,3).*(ones(n_seeds*20,1)*Vbr.dim(1:3));
            br_rand = spm_sample_vol(Vbr,loc_rand(:,1),loc_rand(:,2),loc_rand(:,3),1);
            list = find(br_rand>.9999);
            loc_rand = Vbr.mat*[loc_rand(list,:)' ; ones(1,length(list))];
            loc_rand = loc_rand(1:3,:);
        end
        
        % Remove seeds that are to close (less than 50mm apart) too each other
        if n_seeds>1
            flag=0;
            while flag==0
                dist_loc_r = sqrt(sum((diff(loc_rand').^2)'));
                rm_loc = find(dist_loc_r<50);
                if isempty(rm_loc)
                    flag=1;
                else
                    loc_rand(:,rm_loc) = [];
                end;
            end
            fit_opt.set_loc_o = loc_rand(:,(1:n_seeds+n_dip-1));
        end
    end
    % Put the seeds in spherical coordinates.
    [fit_opt.set_loc_o] = spm_eeg_inv_vert2Rsph(1,model.spheres,fit_opt.set_loc_o);
    
    
    % Prepare head model structure to speed up things.
    %-------------------------------------------------
    MODEL = model;
    MODEL.sigma = model.param.sigma;
    try
       rmfield(MODEL,'flags');
       rmfield(MODEL,'param');
    end
    n_surf = length(model.head);
    MODEL.n_dip = n_dip;
    MODEL.Ntb = size(V,2);
    
    pow2EEG = sum(V_EEG.^2);
    [kk,MODEL.Mtb] =  max(pow2EEG);
    if fit_opt.or_opt==2
        MODEL.orW = kron(sqrt(pow2EEG'),eye(n_dip));
    end
    
    
	% Prepare result files.
	%----------------------
	sdip.n_seeds = n_seeds;
	sdip.n_dip   = n_dip;
	sdip.Sloc    = cell(1,n_seeds);
	sdip.loc     = cell(1,n_seeds);
	sdip.L       = cell(1,n_seeds);
	sdip.j       = cell(1,n_seeds);
	sdip.res     = zeros(1,n_seeds);
	sdip.cost    = zeros(1,n_seeds);
	sdip.M       = Vbr.mat;
	sdip.Mtb     = MODEL.Mtb;
	sdip.fit_opt = fit_opt;
end

% Optimization bit
%-----------------
if ~fit_opt.q_fxd_loc && ~fit_opt.q_fxd_or % Both localisation and orientation have to be optimized
    sdip.exitflag = zeros(1,n_seeds);
	opts = optimset('MaxIter',1000*n_dip,'MaxFunEvals',5000*n_dip);
	% opts = optimset('MaxIter',2000,'MaxFunEvals',10000,'Display','iter','Diagnostics','on')
	for ii=1:n_seeds
        if disp_f, fprintf('\nSeed # %i\n',ii); end
        [sdip.Sloc{ii},fval,sdip.exitflag(ii),outp] = fminsearch('spm_eeg_inv_ecd_costSd',...
                fit_opt.set_loc_o(:,ii+(1:n_dip)-1),opts) ;
        [sdip.cost(ii),sdip.L{ii},sdip.j{ii},sdip.res(ii)] = spm_eeg_inv_ecd_costSd(sdip.Sloc{ii}) ;
        if n_dip==1 && disp_f
            fprintf('\t Sphere location: %4.2f %4.2f %4.2f , source str: %4.2f \n', ...
                [sdip.Sloc{ii}' norm(sdip.j{ii}(:,sdip.Mtb))])
            fprintf('\t Source or: %4.2f %4.2f %4.2f , residual: %4.2f \n', ...
                [sdip.j{ii}(:,sdip.Mtb)'/norm(sdip.j{ii}(:,sdip.Mtb)) sdip.res(ii)])
        end
	end
elseif fit_opt.q_fxd_loc && ~fit_opt.q_fxd_or % Localisation is fixed but not the orientation
    sdip.Sloc{1} = fit_opt.set_loc_o;
    [sdip.cost(1),sdip.L{1},sdip.j{1},sdip.res(1)] = spm_eeg_inv_ecd_costSd(sdip.Sloc{1}) ;
else % Both localisation and orientation are fixed
    sdip.Sloc{1} = fit_opt.set_loc_o;
    [sdip.cost(1),sdip.L{1},sdip.j{1},sdip.res(1)] = spm_eeg_inv_ecd_costSd(sdip.Sloc{1}) ;
end

sdip.rres = sdip.res/norm(V_EEG,'fro'); 
sdip.varexpl = 1-sdip.rres.^2; 

% Turn locations back into real world coordinates
if disp_f
    fprintf('\n\n')
end
for ii=1:n_seeds
    [sdip.loc{ii}] = spm_eeg_inv_vert2Rsph(2,model.spheres,sdip.Sloc{ii});
    if disp_f
        if n_dip==1
            fprintf('\t Real location: %4.2f %4.2f %4.2f \n', sdip.loc{ii}')
        end
        fprintf('\t Relative residual: %4.2f  , Variance explained: %4.2f \n', ...
            [sdip.rres(ii) sdip.varexpl(ii)])
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [c19,L19,j19,res19] = cost_stsm(l19);
% sdip_e19.L{1}


