function D = spm_eeg_inv_vbecd_gui(D,val)
% GUI function for the VB-ECD inversion
% - load the necessary data, if not provided
% - fill in all the necessary bits for the VB-ECD inversion routine,
%   especially the priors
% - launch the VB_ECD routine, aka. spm_eeg_inv_vbecd
% - displays the results.
% See the comments in the function and the article for further details
% 
% Reference:
% Kiebel et al., Variational Bayesian inversion of the equivalent current
% dipole model in EEG/MEG., NeuroImage, 39:728-741, 2008
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips
% $Id: spm_eeg_inv_vbecd_gui.m 3237 2009-06-30 15:20:04Z gareth $

%%
% Load data, if necessary
%==========
if nargin<1
    D = spm_eeg_load;
end

%%
% Check if the forward model was prepared & handle the other info bits
%========================================
if ~isfield(D,'inv')
    error('Data must have been prepared for inversion procedure...')
end
if nargin==2
    % check index provided
    if val>length(D.inv)
        val = length(D.inv);
        D.val = val;
    end
else
    if isfield(D,'val')
        val = D.val;
    else
        % use last one
        val = length(D.inv);
        D.val = val;
    end
end

% Use val to define which is the "current" inv{} to use
% If no inverse solution already calculated (field 'inverse' doesn't exist) 
% use that inv{}. Otherwise create a new one by copying the previous 
% inv{} structure
if isfield(D.inv{val},'inverse')
    % create an extra inv{}
    Ninv = length(D.inv);
    D.inv{Ninv+1} = D.inv{val};
    if isfield(D.inv{Ninv+1},'contrast')
        % no contrast field used here !
        D.inv{Ninv+1} = rmfield(D.inv{Ninv+1},'contrast');
    end
    val = Ninv+1;
    D.val = val;
end

if ~isfield(D.inv{val}, 'date')
    % Set time , date, comments & modality
    clck = fix(clock);
    if clck(5) < 10
        clck = [num2str(clck(4)) ':0' num2str(clck(5))];
    else
        clck = [num2str(clck(4)) ':' num2str(clck(5))];
    end
    D.inv{val}.date = strvcat(date,clck); %#ok<VCAT>
end

if ~isfield(D.inv{val}, 'comment'), 
   D.inv{val}.comment = {spm_input('Comment/Label for this analysis', '+1', 's')};
end

D.inv{val}.method = 'vbecd';

%% Struct that collects the inputs for vbecd code
P = [];

%P.modality = 'EEG';

% Uncomment the line below to try other modalities
P.modality = spm_eeg_modality_ui(D, 1, 1);

if isfield(D.inv{val}, 'forward') && isfield(D.inv{val}, 'datareg')
    for m = 1:numel(D.inv{val}.forward)
        if strncmp(P.modality, D.inv{val}.forward(m).modality, 3)
            P.forward.vol  = D.inv{val}.forward(m).vol;
            if ischar(P.forward.vol)
                P.forward.vol = fileio_read_vol(P.forward.vol);
            end
            P.forward.sens = D.inv{val}.datareg(m).sensors;
            % Channels to use
            P.Ic = setdiff(meegchannels(D, P.modality), badchannels(D));
            
            
            M1 = D.inv{val}.datareg.toMNI;
            if ~isequal(P.modality,'EEG')
                [U, L, V] = svd(M1(1:3, 1:3));
                M1(1:3,1:3) =U*V';
            end
%             disp('Undoing transformation to Tal space !');
%             M1=eye(4)
%             
            P.forward.sens = forwinv_transform_sens(M1, P.forward.sens);
            P.forward.vol = forwinv_transform_vol(M1, P.forward.vol);
            
        end
    end
end

if isempty(P.Ic)
    error(['The specified modality (' P.modality ') is missing from file ' D.fname]);
else
    P.channels = D.chanlabels(P.Ic);
end

 
[P.forward.vol, P.forward.sens] =  forwinv_prepare_vol_sens( ...
    P.forward.vol, P.forward.sens, 'channel', P.channels);

if ~isfield(P.forward.sens,'prj')
    P.forward.sens.prj = D.coor2D(P.Ic);
end


%% 
% Deal with data
%===============

% time bin or time window
msg_tb = ['time_bin or time_win [',num2str(round(min(D.time)*1e3)), ...
            ' ',num2str(round(max(D.time)*1e3)),'] ms'];
ask_tb = 1;
while ask_tb
    tb = spm_input(msg_tb,1,'r');   % ! in msec
    if length(tb)==1
        if tb>=min(D.time([], 'ms')) && tb<=max(D.time([], 'ms'))
            ask_tb = 0;
        end
    elseif length(tb)==2
        if all(tb>=floor(min(D.time([], 'ms')))) && all(tb<=ceil(max(D.time([], 'ms')))) && tb(1)<=tb(2)
            ask_tb = 0;
        end
    end
end
if length(tb)==1
    [kk,ltb] = min(abs(D.time([], 'ms')-tb)); % round to nearest time bin
else
    [kk,ltb(1)] = min(abs(D.time([], 'ms')-tb(1)));  % round to nearest time bin
    [kk,ltb(2)] = min(abs(D.time([], 'ms')-tb(2)));
    ltb = ltb(1):ltb(2); % list of time bins 'tb' to use
end

%% get a baseline period- used to get precision
msg_tb = ['Baseline time_win [',num2str(round(min(D.time)*1e3)), ...
            ' ',num2str(round(max(D.time)*1e3)),'] ms'];
ask_tb = 1;
while ask_tb
    btb = spm_input(msg_tb,1,'r');   % ! in msec
    if length(btb)==1
        if btb>=min(D.time([], 'ms')) && btb<=max(D.time([], 'ms'))
            ask_tb = 0;
        end
    elseif length(btb)==2
        if all(btb>=floor(min(D.time([], 'ms')))) && all(btb<=ceil(max(D.time([], 'ms')))) && btb(1)<=btb(2)
            ask_tb = 0;
        end
    end
end
if length(btb)==1
    [kk,bltb] = min(abs(D.time([], 'ms')-btb)); % round to nearest time bin
else
    [kk,bltb(1)] = min(abs(D.time([], 'ms')-btb(1)));  % round to nearest time bin
    [kk,bltb(2)] = min(abs(D.time([], 'ms')-btb(2)));
    bltb = bltb(1):bltb(2); % list of time bins 'tb' to use
end

% trial type
if D.ntrials>1
    msg_tr = ['Trial type number [1 ',num2str(D.ntrials),']'];
    ltr = spm_input(msg_tr,2,'i',num2str(1:D.ntrials));
    tr_q = 1;
else
    tr_q = 0;
    ltr = 1;
end

% data, averaged over time window considered


%% convert to VOLTS
if strcmp(upper(P.modality),'EEG'),
   
    
    eegunits = unique(D.units(D.meegchannels('EEG')));
    Neegchans=numel(D.units(D.meegchannels('EEG')));
if ~strcmp(eegunits,'V'),
    warning('units unspecified');
    if mean(std(D(P.Ic,ltb,ltr)))>1e-2,
        guess_units='uV';
        else
        guess_units='V';
        end;
     fprintf('No units specified but rms of data is %3.2f units\n',mean(std(D(P.Ic,ltb,ltr))))
     msg_str=sprintf('Units of EEG are %s ?',guess_units);
    ans=spm_input(msg_str,1,'s','yes');   
    if ~strcmp(ans,'yes')
        error('stop for now');
        end; % if strcmp
    D = units(D, 1:Neegchans, guess_units)
    end; %if no units
    
    eegunits = unique(D.units(D.meegchannels('EEG')));
    EEGscale=-1;
    if strcmp(eegunits,'V'),
            EEGscale=1;
    end; % if
    if strcmp(eegunits,'uV'),
            EEGscale=1e-6;
    end % if
    if EEGscale==-1,
        error('unknown eeg units');
    end; % if
end; % if eeg data

EEGscale
dat_y = squeeze(mean(D(P.Ic,ltb,ltr)*EEGscale,2));
base_dat_stdy = squeeze(std(D(P.Ic,bltb,ltr)*EEGscale,0,2)); %% baseline data- standard deviation

%%
% Other bits of the P structure, apart for priors and #dipoles
%==============================

P.ltr          = ltr;
P.Nc           = length(P.Ic);
P.Niter        = 200;           % \_ Using SK default values here...
P.threshold_dF = 1e-2;  %1e-4   % /

%%  
% Deal with dipoles number and priors
%====================================
dip_q = 0; % number of dipole 'elements' added (single or pair)
dip_c = 0; % total number of dipoles in the model
adding_dips = 1;
clear dip_pr

% These are defaults value picked from SK's example data set ran for the
% paper, values should proabably be tweaked a bit...
%def_ab_noninfo = [1e-3 1e-12];
%def_ab_info    = [3 1e-12];

%% new default priors
%% b20/a20 should be the expect error variance in source moment (in nAM ?). 
%% a,b determine distribution of hyper-priors (which determine paramter
%% variances)
%% b/a gives the expected mean of this dist'n, i.e. the variance of the
%% prior
orient_noninfo=[1.5 15]; %% expected variance is given by 15/1.5
orient_info=[1.5 0.15];

pos_noninfo=[1.5 500];  %% approx variability of sqrt(500) mm
pos_info=[1.5 100];  %% approx variability of sqrt(100) mm

%% b30/a30 should be the expected error variance in location parameters (in
%% mm) 


corr = .999;

while adding_dips
    if dip_q>0, 
        msg_dip =['Add dipoles to ',num2str(dip_c),' or stop?'];
        dip_ch = 'Single|Pair|Stop';
        dip_val = [1,2,0];
        def_opt=3;
    else
        msg_dip =['Add dipoles to model'];
        def_opt=1;
        dip_ch = 'Single|Pair';
        dip_val = [1,2];
    end
    a_dip = spm_input(msg_dip,2+tr_q+dip_q,'b',dip_ch,dip_val,def_opt);
    if a_dip == 0
        adding_dips = 0;
    elseif a_dip == 1
    % add a single dipole to the model
        dip_q = dip_q+1;
        dip_pr(dip_q) = struct( 'a_dip',a_dip, ...
            'mu_w0',[],'mu_s0',[],'S_s0',eye(3),'S_w0',eye(3), ...
            'ab20',[],'ab30',[],'Tw',eye(3),'Ts',eye(3));
        % Location prior
        spr_q = spm_input('Location prior ?',1+tr_q+dip_q+1,'b', ...
                    'Informative|Non-info',[1,0],2);
        if spr_q
            % informative location prior
            str = 'Location prior';
            while 1
                s0 = spm_input(str, 1+tr_q+dip_q+2,'e',[0 0 0])';
                outside = ~forwinv_inside_vol(s0',P.forward.vol);
                if all(~outside), break, end
                str = 'Prior location must be inside head';
            end
            dip_pr(dip_q).mu_s0 = s0;
            dip_pr(dip_q).ab30 = pos_info;
        else
            % no location  prior
            dip_pr(dip_q).mu_s0 = zeros(3,1);
            dip_pr(dip_q).ab30 = pos_noninfo;
        end
        % Moment prior
        wpr_q = spm_input('Moment prior ?',1+tr_q+dip_q+spr_q+2,'b', ...
                    'Informative|Non-info',[1,0],2);
        if wpr_q
            % informative moment prior
            dip_pr(dip_q).mu_w0 = spm_input('Moment prior', ...
                                        1+tr_q+dip_q+spr_q+3,'e',[0 0 0])';
            dip_pr(dip_q).ab20 = orient_info;
        else
            % no location  prior
            dip_pr(dip_q).mu_w0 = zeros(3,1);
            dip_pr(dip_q).ab20 = orient_noninfo;
        end
        dip_c = dip_c+1;
    else
    % add a pair of symmetric dipoles to the model
        dip_q = dip_q+1;
        dip_pr(dip_q) = struct( 'a_dip',a_dip, ...
            'mu_w0',[],'mu_s0',[],'S_s0',eye(6),'S_w0',eye(6), ...
            'ab20',[],'ab30',[],'Tw',eye(6),'Ts',eye(6));
        % Location prior
        spr_q = spm_input('Location prior ?',1+tr_q+dip_q+1,'b', ...
                    'Informative|Non-info',[1,0],2);
        if spr_q
            % informative location prior
            str = 'Location prior (right only)';
            while 1
                tmp = spm_input(str, 1+tr_q+dip_q+2,'e',[0 0 0])';
                outside = ~forwinv_inside_vol(tmp',P.forward.vol);
                if all(~outside), break, end
                str = 'Prior location must be inside head';
            end
            tmp = [tmp ; tmp] ; tmp(4) = -tmp(4);
            dip_pr(dip_q).mu_s0 = tmp ;
            dip_pr(dip_q).ab30 = pos_info;
        else
            % no location  prior
            dip_pr(dip_q).mu_s0 = zeros(6,1);
            dip_pr(dip_q).ab30 = pos_noninfo;
        end
        % Moment prior
        wpr_q = spm_input('Moment prior ?',1+tr_q+dip_q+spr_q+2,'b', ...
                                           'Informative|Non-info',[1,0],2);
        if wpr_q
            % informative moment prior
            tmp = spm_input('Moment prior (right only)', ...
                                        1+tr_q+dip_q+spr_q+3,'e',[1 1 1])';
            tmp = [tmp ; tmp] ; tmp(4) = -tmp(4);
            dip_pr(dip_q).mu_w0 = tmp;
            dip_pr(dip_q).ab20 = orient_info;
        else
            % no location  prior
            dip_pr(dip_q).mu_w0 = zeros(6,1);
            dip_pr(dip_q).ab20 = orient_noninfo;
        end
        % Symmetry priors, soft or hard for both location and moment.
        mpr_q = spm_input('Symmetry prior ?', ...
                       1+tr_q+dip_q+spr_q+wpr_q+3,'b','Soft|Hard',[1,0],1);
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
                'a20',[],'b20',[], ...
                'a30',[],'b30',[]);
% PROBLEM:
% How to ensure different precision level, i.e. informative vs
% non-informative, for different priors !!!
% Trying to impose that with an a priori scaling of the prior covariance 
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
    for ii=1:dip_c
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

% Initialise inverse field
inverse = struct( ...
    'F',[], ... % free energy
    'pst',D.time, ... % all time points in data epoch
    'tb',tb, ... % time window/bin used
    'ltb',ltb, ... % list of time points used
    'ltr',ltr, ... % list of trial types used
    'n_seeds',length(ltr), ... % using this field for multiple reconstruction
    'n_dip',dip_c, ... % number of dipoles used
    'loc',[], ... % loc of dip (3 x n_dip)
    'j',[], ... % dipole(s) orient/ampl, in 1 column
    'cov_loc',[], ... % cov matrix of source location
    'cov_j',[], ...   % cov matrix of source orient/ampl
    'Mtb',1, ... % ind of max EEG power in time series, 1 as only 1 tb.
    'exitflag',[], ... % Converged (1) or not (0)
    'P',[]);             % save all kaboodle too.

for ii=1:length(ltr)
    P.y = dat_y(:,ii);
    P.ii = ii;
    
   




%   disp('FIXING DIP TO GET SCALING');
%   mu_s=[-30 20 20]';
%   mu_w0=[ 1 1 0]';
%   P.dv=10^-2;
%   
%   [gmn, gm, dgm] = spm_eeg_inv_vbecd_getLF(mu_s, P.forward.sens, P.forward.vol,...
%        P.dv.*ones(1, length(mu_w0))); %% order of 10-4 for eeg
%  
%signalmag=max(std(dat_y));
%disp('Guessing at signal magnitude and SNR');
%signalmag=5e-4;

chanSNR=(dat_y./(base_dat_stdy)).^2;
powerSNR=mean(chanSNR);
disp(sprintf('Estimating SNR (power) to be %3.2f',powerSNR));
errorvar=(mean(base_dat_stdy)).^2; %% variance of error 
alla10 = numel(P.y)/2; 
allb10 = alla10*errorvar/2; % data
P.priors.a10=alla10;
P.priors.b10=allb10;

testing=0;
if testing,
   
%% b10/a10 should be the expected error variance at the sensors -ideally
%% should come from the standard error of the average trace

%% b20/a20 should be the expect error variance in source moment (in nAM ?). 
alla20=1.5; %% half number of orientation params for now
allb20=alla20*10; %% 

%% b30/a30 should be the expected error variance in location parameters (in
%% mm) 
alla30=1.5; %% half num of position params
allb30=100;

P.priors.a20=alla20;
%% prior info on b20
P.priors.b20=allb20;
P.priors.b30=allb30;
P.priors.a30=alla30;
end; % if 

    fprintf('\nLocalising source nr %d.\n',ii)

    P = spm_eeg_inv_vbecd(P);

  
    ssres=sum(P.post.residuals.^2)./length(P.y);
    disp(sprintf('prior estimate var=%3.2e, post=  %3.2e, actual resid=%3.2e',P.priors.b10/P.priors.a10,P.post.b1/P.post.a1,ssres));
    disp(sprintf('prior estimate moment var=%3.2f , post=  %3.2f',P.priors.b20/P.priors.a20,P.post.b2/P.post.a2));
    disp(sprintf('prior estimate pos var=%3.2f (=%3.2fmm per dim), post= %3.2f',P.priors.b30/P.priors.a30,power(P.priors.b30/P.priors.a30,1/2),P.post.b3/P.post.a3));
    P.post.mu_s'
    P.post.mu_w'
 
    % Get the results out.
    inverse.pst = tb*1e3;
    inverse.F(ii) = P.F; % free energy
    inverse.loc{ii} = reshape(P.post.mu_s,3,dip_c); % loc of dip (3 x n_dip)
    inverse.j{ii} = P.post.mu_w; % dipole(s) orient/ampl, in 1 column
    inverse.cov_loc{ii} = P.post.S_s; % cov matrix of source location
    inverse.cov_j{ii} = P.post.S_w; % cov matrix of source orient/ampl
    inverse.exitflag(ii) = P.ok; % Converged (1) or not (0)
    inverse.P{ii} = P; % save all kaboodle too.
end
D.inv{val}.inverse = inverse;

%%
% Save results and display
%-------------------------
save(D)
P.handles.hfig  = spm_figure('GetWin','Graphics');
close(P.handles.hfig)
spm_eeg_inv_vbecd_disp('init',D);

return
