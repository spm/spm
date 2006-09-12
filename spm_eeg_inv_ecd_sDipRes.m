function [result,Pres] = spm_eeg_inv_ecd_sDipRes(sdip)

%____________________________________________________________________________
%
% FORMAT [result,Pres] = spm_eeg_inv_ecd_sDipRes(sdip)
%
% This function summarizes the result of the dipole fitting routines.
% Dipoles are grouped by clusters (or group) according to their location.
%
% The grouping critiria is based on the relative residual error (rres):
%   if many sets have very similar 'rres', then their dipoles location
%   is also similar.
% Then these locations are averaged to give the "center" of the cluster,
% the other infos (j, res, rres, cost, etc) are also averaged.
% Finally clusters are ordered according to the number of sets they summarize.
%
% To display the results, I can use:
%       spm_eeg_inv_ecd_DrawDip('init',result)
%____________________________________________________________________________
%
% 'result' is a structure with following fields
% - lrres: list of results, ordered according to their 'rres'.
% - n_gr: # of source sets per group.
% - N_gr: # of groups comprising similar fitted sources
% - perm: permutation of sources within a set, cell{1 x N_gr}(n_gr x 2)
% - dmin: distance between the dipoles of 1st set and other sets inthe group,
%         cell{1 x N_gr}(n_gr x 2)
% - j: "Mean" amplitude of dipoles in the group, cell{1 x N_gr}(2*n_dip x Ntimebins)
% - j3: "Mean" amplitude of dipoles in the group, cell{1 x N_gr}(3 x Ntimebins x n_dip)
% - or: "Mean" orientation of dipoles in the group, cell{1 x N_gr}(3 x n_dip)
% - loc: "Mean" location of dipoles in the group, cell{1 x N_gr}(3 x n_dip)
% - res: "Mean" residuals of dipole fit, vect(1 x N_gr)
% - rres: "Mean" relative residuals of dipole fit, vect(1 x N_gr)
% - cost: "Mean" cost of dipole fit, vect(1 x N_gr)
% - varexpl: variance explained byt fitted dipole(s)
% - n_seeds: # of groups left, similar to # of seeding sets
% - n_dip: # of dipoles per fit
% - exitflag: how the fitting process converged
% - Mtb: index of maximum power in EEG time series used
% - M: 4x4 affine transformation matrix mapping from brain image voxel coordinates
%            to real world coordinates
% Fields 'n_seeds', 'n_dip', 'exitflag', 'Mtb', 'M' are added so that I can display
% the 'result' with draw_sdip
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips,
% $Id$

if nargin < 1
    Pdip = spm_select(1,'^S.*dip.*\.mat$','Select dipole file');
    load(Pdip);
    flag_file = 1;
else
    flag_file = 0;
end


% Sort solutions according to their relative residuals
%=====================================================
[srres,lrres] = sort(sdip.rres);
rdsrres = diff(srres/min(srres)*100);
l_gr   = find(rdsrres>=.1);
ll_gr  = [0 l_gr];
n_gr   = diff([ll_gr sdip.n_seeds]);
N_gr   = length(n_gr);

[sn_gr,perm_gr] = sort(-n_gr);
sn_gr = -sn_gr;
kk    = [];
for ii=1:N_gr
    l_g = ll_gr(perm_gr(ii))+(1:n_gr(perm_gr(ii)));
    kk = [kk lrres(l_g)];
end
lrres = kk;
n_gr  = sn_gr;
ll_gr = [0 cumsum(n_gr)];

% Prepare result structure
%=========================
result.lrres = lrres;
result.n_gr = n_gr;
result.N_gr = N_gr;
result.perm = cell(1,N_gr); % order of the dipoles compared to the 1st in the group
result.dmin = cell(1,N_gr); % Distance beween the 1st dipoles and the ones from other groups
result.j = cell(1,N_gr);    % "Mean" amplitude of dipoles in the group, as 3n_dip x Ntb 
result.j3 = cell(1,N_gr);   % "Mean" amplitude of dipoles in the group, as 3 x Ntb x n_dip
result.or = cell(1,N_gr);   % "Mean" orientation of dipoles in the group
result.loc = cell(1,N_gr);  % "Mean" location of dipoles in the group
result.L = cell(1,N_gr);    % "Mean" leadfield of dipoles in the group
result.res = zeros(1,N_gr);
result.rres = zeros(1,N_gr);
result.cost = zeros(1,N_gr);
result.varexpl = zeros(1,N_gr);

Ntb = size(sdip.j{1},2);
unit = ones(1,sdip.n_dip);
shift = [-2 -1 0 -2 -1 0];
% Deal with the various 'groups' of solutions.
for ii=1:N_gr
    l_g = ll_gr(ii)+(1:n_gr(ii));
    if n_gr(ii)==1
        % Lonely solution
        result.j3{ii} = zeros(3,Ntb,sdip.n_dip);
        for jj=1:sdip.n_dip
            result.j3{ii}(:,:,jj) = sdip.j{lrres(l_g(1))}(3*jj+(-2:0),:);
        end
        result.j{ii} = sdip.j{lrres(l_g(1))};
        result.loc{ii} = sdip.loc{lrres(l_g(1))};
        result.L{ii} = sdip.L{lrres(l_g(1))};
        result.res(ii) = sdip.res(lrres(l_g(1)));
        result.rres(ii) = sdip.rres(lrres(l_g(1)));
        result.cost(ii) = sdip.cost(lrres(l_g(1)));
   else
        % Determine the permutation such that in a group dipoles correspond to each other
        % The 1st dipole set of the group is used as reference.
        result.perm{ii} = zeros(n_gr(ii),sdip.n_dip);
        result.dmin{ii} = zeros(n_gr(ii),sdip.n_dip);
        result.perm{ii}(1,:) = 1:sdip.n_dip ; 
        for jj=2:n_gr(ii)
            d_loc = zeros(sdip.n_dip);
            for kk=1:sdip.n_dip
                d_loc(:,kk) = sqrt( sum( ...
                    (sdip.loc{lrres(l_g(1))}(:,kk)*unit - sdip.loc{lrres(l_g(jj))}).^2 )' );
            end
            for kk=1:sdip.n_dip     
                [vmin,pmin] = min(d_loc(:));
                n = floor((pmin-1)/sdip.n_dip)+1;
                m = mod(pmin-1,sdip.n_dip)+1;
                result.perm{ii}(jj,n) = m ;
                result.dmin{ii}(jj,n) = vmin ;
                d_loc(m,:) = 1e6; d_loc(:,n) = 1e6;
            end
        end
        % Mean loc and amplitudes
        result.loc{ii} = zeros(3,sdip.n_dip);
        result.L{ii} = zeros(size(sdip.L{1}));
        j_tmp = zeros(size(sdip.j{ii}));
        pp = kron(ones(1,sdip.n_dip),-2:0);
        for jj=1:n_gr(ii)
            result.loc{ii} = result.loc{ii} + sdip.loc{lrres(l_g(jj))}(:,result.perm{ii}(jj,:));
            result.L{ii} = result.L{ii} + ...
                sdip.L{lrres(l_g(jj))}(:,kron(3*result.perm{ii}(jj,:),[1 1 1])+pp);
            j_tmp = j_tmp + ...
                sdip.j{lrres(l_g(jj))}(kron(result.perm{ii}(jj,:)*3,[1 1 1])+pp,:);
        end
        result.loc{ii} = result.loc{ii}/n_gr(ii);
        result.L{ii} = result.L{ii}/n_gr(ii);
        result.j3{ii} = zeros(3,Ntb,sdip.n_dip);
        for jj=1:sdip.n_dip
            result.j3{ii}(:,:,jj) = j_tmp(3*jj+(-2:0),:)/n_gr(ii);
        end
        result.j{ii} = j_tmp/n_gr(ii);
        result.res(ii) = mean(sdip.res(lrres(l_g)));
        result.rres(ii) = mean(sdip.rres(lrres(l_g)));
        result.cost(ii) = mean(sdip.cost(lrres(l_g)));
    end
    or_tmp = result.j3{ii}./repmat(sqrt(sum(result.j3{ii}.^2)),[3,1,1]);
    if size(or_tmp,2)==1
       result.or{ii} = or_tmp;
   else
       if ~isfield(sdip,'fit_opt')
          result.or{ii} = or_tmp;
       elseif sdip.fit_opt.or_opt==2 | sdip.fit_opt.or_opt==3
          result.or{ii} = squeeze(or_tmp(:,1,:));
       elseif sdip.fit_opt.or_opt==4 & sdip.fit_opt.q_fxd_or
          result.or{ii} = squeeze(or_tmp(:,1,:));
       else
          result.or{ii} = or_tmp;
       end
  end
end
result.varexpl = 1-result.rres.^2;


% Adapt things so that I can display the 'result' with draw_sdip
%===============================================================
result.n_seeds = N_gr;
result.n_dip = sdip.n_dip;
result.exitflag = ones(1,N_gr);
result.Mtb = sdip.Mtb;
result.M = sdip.M;
result.fit_opt = sdip.fit_opt;
if isfield(sdip,'Pdata')
    result.Pdata = sdip.Pdata;
elseif isfield(sdip,'data')
    result.Pdata = sdip.data;
end

% Save things if necessary & display
if flag_file
    resdip = result;
    Pres = [spm_str_manip(Pdip,'H'),filesep,'res_',spm_str_manip(Pdip,'t')];
    save(Pres,'resdip')
else
    Pres = '';
end

return

