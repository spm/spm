% spm_eeg_inv_vbecd_gui
%
% Script function to 
% - load data
% - fill in all the necessary bits for the VB-ECD inversion routine
% - launch the VB_ECD routine
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips
% $Id: spm_eeg_inv_vbecd_gui.m 1932 2008-07-18 17:12:02Z christophe $

%%
% Load data
%==========
D = spm_eeg_load;

%%
% Check if the forward model was prepared
%========================================
if ~isfield(D,'inv')
    error('Data must have been prepared for inversion procedure...')
else
    if isfield(D,'val')
        val = D.val;
    else
        % use last one
        val = length(D.inv);
        D.val = val;
    end
end

if isfield(D.inv{val},'forward')
    if isfield(D.inv{val}.forward,'vol')
        P.forward = struct('vol',D.inv{val}.forward.vol, ...
                            'sens',[]);
    else
        error('Forward model needs to be ready in FT format.!')
    end
    if isfield(D.inv{val}.forward,'sens')
        P.forward.sens = D.inv{val}.forward.sens;
    else
        P.forward.sens = D.sensors('eeg');
    end
else
    error('Forward model needs to be ready in FT format.!')
end


%% 
% Deal with data
%===============
P.Bad = D.badchannels;

% time bin or time window
msg_tb = ['t_b or t_w [',num2str(min(D.time)*1000), ...
            ' ',num2str(max(D.time)*1000),'] ms'];
ask_tb = 1
while ask_tb
    tb = spm_input(msg_tb,1,'r')/1000;
    if length(tb)==1
        if tb>=min(D.time) && tb<=max(D.time)
            ask_tb = 0;
        end
    elseif length(tb)==2
        if all(tb>=min(D.time)) && all(tb<=max(D.time)) && tb(1)<=tb(2)
            ask_tb = 0;
        end
    end
end
if length(tb)==1
    [kk,ltb] = min(abs(D.time-tb)); % round to nearest time bin
else
    [kk,ltb(1)] = min(abs(D.time-tb(1)));  % round to nearest time bin
    [kk,ltb(2)] = min(abs(D.time-tb(2)));
    ltb = ltb(1):ltb(2); % list of time bins 'tb' to use
end

% trial type
if D.ntrials>1
    msg_tr = ['Trial type number [1 ',num2str(D.ntrials),']'];
    ltr = spm_input(msg_tr,2,'i','1');
else
    ltr=1;
end

% data, averaged over time window considered
P.y = mean(squeeze(D(meegchannels(D),ltb,ltr)),2);

%%
% Other bits of the P structure, apart for priors and #dipoles
%==============================

P.Nc           = length(meegchannels(D))-length(P.Bad);
P.Niter        = 200;  % \_ Using SK default values here...
P.threshold_dF = 1e-4; % /

%%  
% Deal with dipoles number and priors
%====================================
dip_q = 1;
dip_c = 0;
adding_dips = 1;
clear dip_pr

% These are defaults value picked from SK's example data set ran for the
% paper, values should proabably be tweaked a bit...
def_ab_noninfo = [1e-3 1e-12];
def_ab_info    = [3 1e-12];
corr = .999;

while adding_dips
    msg_dip =['Add dipoles to ',num2str(dip_c),'  or stop?'];
    a_dip = spm_input(msg_dip,2+dip_q,'b','Single|Pair|Stop',[1,2,0],1);
    if a_dip == 0
        adding_dips = 0;
    elseif a_dip == 1
        % add a single dipole to the model
        dip_pr(dip_q) = struct( 'a_dip',a_dip, ...
            'mu_w0',[],'mu_s0',[],'S_s0',eye(3),'S_w0',eye(3), ...
            'ab20',[],'ab30',[],'Tw',eye(3),'Ts',eye(3));
        % Location prior
        spr_q = spm_input('Location prior ?',2+dip_q,'b', ...
                    'Informative|Non-info',[1,0],2);
        if spr_q
            % informative location prior
            dip_pr(dip_q).mu_s0 = spm_input('Location prior',2+dip_q+1, ...
                                    'e',[0 0 0])';
            dip_pr(dip_q).ab30 = def_ab_info;
        else
            % no location  prior
            dip_pr(dip_q).mu_s0 = zeros(3,1);
            dip_pr(dip_q).ab30 = def_ab_noninfo;
        end
        % Moment prior
        wpr_q = spm_input('Moment prior ?',2+dip_q+1+spr_q,'b', ...
                    'Informative|Non-info',[1,0],2);
        if wpr_q
            % informative moment prior
            dip_pr(dip_q).mu_w0 = spm_input('Moment prior',2+dip_q+2+spr_q, ...
                                    'e',[0 0 0])';
            dip_pr(dip_q).ab20 = def_ab_info;
        else
            % no location  prior
            dip_pr(dip_q).mu_w0 = zeros(3,1);
            dip_pr(dip_q).ab20 = def_ab_noninfo;
        end
        dip_c = dip_c+1;
        dip_q = dip_q+1;
    else
        % add a pair of symmetric dipoles to the model
        dip_pr(dip_q) = struct( 'a_dip',a_dip, ...
            'mu_w0',[],'mu_s0',[],'S_s0',eye(6),'S_w0',eye(6), ...
            'ab20',[],'ab30',[],'Tw',eye(6),'Ts',eye(6));
        % Location prior
        spr_q = spm_input('Location prior ?',2+dip_q,'b', ...
                    'Informative|Non-info',[1,0],2);
        if spr_q
            % informative location prior
            tmp = spm_input('Location prior (right only)',...
                                    2+dip_q+1,'e',[10 0 0])';
            tmp = [tmp ; tmp] ; tmp(4) = -tmp(4);
            dip_pr(dip_q).mu_s0 = tmp ;
            dip_pr(dip_q).ab30 = def_ab_info;
        else
            % no location  prior
            dip_pr(dip_q).mu_s0 = zeros(6,1);
            dip_pr(dip_q).ab30 = def_ab_noninfo;
        end
        % Moment prior
        wpr_q = spm_input('Moment prior ?',2+dip_q+spr_q,'b', ...
                    'Informative|Non-info',[1,0],2);
        if wpr_q
            % informative moment prior
            tmp = spm_input('Moment prior (right only)',2+dip_q+2+spr_q, ...
                                    'e',[1 1 1])';
            tmp = [tmp ; tmp] ; tmp(4) = -tmp(4);
            dip_pr(dip_q).mu_w0 = tmp;
            dip_pr(dip_q).ab20 = def_ab_info;
        else
            % no location  prior
            dip_pr(dip_q).mu_w0 = zeros(6,1);
            dip_pr(dip_q).ab20 = def_ab_noninfo;
        end
        % Symmetry priors, soft or hard for both location and moment.
        mpr_q = spm_input('Symmetry prior ?',2+dip_q+3+spr_q+wpr_q,'b', ...
                    'Soft|Hard',[1,0],1);
        if mpr_q
            % Soft prior, i.e. parameter correlation
            tmp = eye(6);
            tmp2 = eye(3); tmp2(1,1) = -1;
            tmp(4:6,1:3) = tmp2*corr; 
            tmp(1:3,4:6) = tmp2*corr; 
            dip_pr(dip_q).S_s0 = tmp;
            dip_pr(dip_q).S_w0 = tmp;
        else
            % hard prior, i.e. parameter reduction
            T = [eye(3) ; eye(3)]; T(4,1) = -1;
            dip_pr(dip_q).Tw = .5*T*T';
            dip_pr(dip_q).Ts = .5*T*T';            
        end
        dip_q = dip_q+1;
        dip_c = dip_c+2;
    end
end

%%
% Get all the priors together and build structure to pass to inv_vbecd !
%============================
priors = struct('mu_w0',cat(1,dip_pr(:).mu_w0), ...
                'mu_s0',cat(1,dip_pr(:).mu_s0), ...
                'iS_w0',[],'iS_s0',[], ...
                'Tw',blkdiag(dip_pr(:).Tw), ...
                'Ts',blkdiag(dip_pr(:).Ts),...
                'a10',1e-3,'b10',1e-12, ...
                'a20',[],'b20',[], ...
                'a30',[],'b30',[]);
% PROBLEM:
% How to ensure different precision level, i.e. informative vs
% non-informative, for different priors !!!
% Trying to impose that with an prior scaling of the prior covariance 
% matrices but I still think this is not the right way to go...

tmp_ab20 = cat(1,dip_pr(:).ab20);
% assume that parameter b20/30 have always the same value...
if length(unique(tmp_ab20(:,1)))==1
    % all model element have same prior precision, easy !
    priors.a20   = unique(tmp_ab20(:,1));
    priors.b20   = unique(tmp_ab20(:,2));   
    priors.iS_w0 = inv(blkdiag(dip_pr(:).S_w0));
else
    % Not equal "informativeness" -> tricky
    % assume 2 levels, and use their ratio to weight variance mtx
    priors.a20   = min(tmp_ab20(:,1));
    priors.b20   = min(tmp_ab20(:,2));   
    for ii=1:dip_q
        tmp_iS{ii} = inv(dip_pr(ii).S_w0*min(tmp_ab20(:,1))/tmp_ab20(ii,1));
    end
    priors.iS_w0 = blkdiag(tmp_iS{:});
end

tmp_ab30 = cat(1,dip_pr(:).ab30);
if length(unique(tmp_ab30(:,1)))==1
    % all model element have same prior precision, easy !
    priors.a30   = unique(tmp_ab30(:,1));
    priors.b30   = unique(tmp_ab30(:,2));   
    priors.iS_s0 = inv(blkdiag(dip_pr(:).S_s0));
else
    % Not equal "informativeness" -> tricky
    % assume 2 levels, and use their ratio to weight variance mtx
    priors.a30   = min(tmp_ab30(:,1));
    priors.b30   = min(tmp_ab30(:,2));   
    for ii=1:dip_q
        tmp_iS{ii} = inv(dip_pr(ii).S_s0*min(tmp_ab30(:,1))/tmp_ab30(ii,1));
    end
    priors.iS_s0 = blkdiag(tmp_iS{:});
end
P.priors = priors;

%%
% Launch inversion !
%===================
P = spm_eeg_inv_vbecd(P)

