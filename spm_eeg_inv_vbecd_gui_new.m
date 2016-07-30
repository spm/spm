function D = spm_eeg_inv_vbecd_gui(D,val)
% GUI function for Bayesian ECD inversion
% - load the necessary data, if not provided
% - fill in all the necessary bits for the VB-ECD inversion routine,
% - launch the B_ECD routine, aka. spm_eeg_inv_vbecd
% - displays the results.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
%
% $Id: spm_eeg_inv_vbecd_gui_new.m 6847 2016-07-30 10:35:32Z karl $

% Load data, if necessary
%==========================================================================
if nargin < 1
    D = spm_eeg_load;
end


% Check if the forward model was prepared & handle the other info bits
%==========================================================================
if ~isfield(D,'inv')
    error('Data must have been prepared for inversion procedure...')
end
if nargin == 2
    % check index provided
    %----------------------------------------------------------------------
    if val>length(D.inv)
        val   = length(D.inv);
        D.val = val;
    end
else
    if isfield(D,'val')
        val = D.val;
    else
        % use last one
        val   = length(D.inv);
        D.val = val;
    end
end

% Use val to define which is the "current" inv{} to use
% If no inverse solution already calculated (field 'inverse' doesn't exist)
% use that inv{}. Otherwise create a new one by copying the previous
%--------------------------------------------------------------------------
if isfield(D.inv{val},'inverse')
    Ninv          = length(D.inv);
    D.inv{Ninv+1} = D.inv{val};
    if isfield(D.inv{Ninv+1},'contrast')
        D.inv{Ninv + 1} = rmfield(D.inv{Ninv + 1},'contrast');
    end
    val   = Ninv + 1;
    D.val = val;
end

% Set time , date, comments & modality
%--------------------------------------------------------------------------
if ~isfield(D.inv{val}, 'date')
    clck = fix(clock);
    if clck(5) < 10
        clck = [num2str(clck(4)) ':0' num2str(clck(5))];
    else
        clck = [num2str(clck(4)) ':' num2str(clck(5))];
    end
    D.inv{val}.date = strvcat(date,clck);
end

if ~isfield(D.inv{val}, 'comment'),
    D.inv{val}.comment = {spm_input('Comment/Label for this analysis', '+1', 's')};
end

D.inv{val}.method = 'vbecd';

% Struct that collects the inputs for vbecd code
%--------------------------------------------------------------------------
P          = [];
P.modality = spm_eeg_modality_ui(D, 1, 1);

for m = 1:numel(D.inv{val}.datareg)
    if strncmp(D.inv{val}.datareg(m).modality, P.modality, 3)
        modind = m;
    end
end

data              = spm_eeg_inv_get_vol_sens(D, val, [], 'inv', P.modality);
P.forward.vol     = data.(P.modality(1:3)).vol;
if ischar(P.forward.vol)
    P.forward.vol = ft_read_vol(P.forward.vol);
end
P.forward.sens    = data.(P.modality(1:3)).sens;
P.forward.siunits = data.siunits;
M1 = data.transforms.toMNI;
[U,L,V] = svd(M1(1:3, 1:3));
orM1(1:3,1:3) = U*V'; % for switching orientation between meg and mni space

P.Ic = indchantype(D, P.modality, 'GOOD');

if isempty(P.Ic)
    error(['The specified modality (' P.modality ') is missing from file ' D.fname]);
else
    P.channels = D.chanlabels(P.Ic);
end

P.forward.chanunits = D.units(P.Ic);
[P.forward.vol, P.forward.sens] =  ft_prepare_vol_sens( ...
    P.forward.vol, P.forward.sens, 'channel', P.channels);
P.forward.sens.prj  = D.coor2D(P.Ic);

% Deal with data
%==========================================================================

% time bin or time window
%--------------------------------------------------------------------------
msg_tb = ['time_bin or average_win [',num2str(round(min(D.time)*1e3)), ...
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
if length(tb) == 1
    [kk,ltb]    = min(abs(D.time([], 'ms')-tb));     % round to nearest time bin
else
    [kk,ltb(1)] = min(abs(D.time([], 'ms')-tb(1)));  % round to nearest time bin
    [kk,ltb(2)] = min(abs(D.time([], 'ms')-tb(2)));
    ltb         = ltb(1):ltb(2);                     % list of time bins 'tb' to use
end

% trial type
%--------------------------------------------------------------------------
if D.ntrials>1
    msg_tr = ['Trial type number [1 ',num2str(D.ntrials),']'];
    ltr    = spm_input(msg_tr,2,'i',num2str(1:D.ntrials));
    tr_q   = 1;
else
    tr_q   = 0;
    ltr    = 1;
end

% data, averaged over time window considered
%--------------------------------------------------------------------------

EEGscale = 1;

% SORT OUT EEG UNITS AND CONVERT VALUES TO VOLTS
%--------------------------------------------------------------------------
if strcmpi(P.modality,'EEG')
    allunits  = strvcat('uV','mV','V');
    allscales = [1, 1e3, 1e6]; %
    EEGscale  = 0;
    eegunits  = unique(D.units(D.indchantype('EEG')));
    Neegchans = numel( D.units(D.indchantype('EEG')));
    for j=1:length(allunits)
        if strcmp(deblank(allunits(j,:)),deblank(eegunits))
            EEGscale=allscales(j);
        end % if
    end % for j
    
    if EEGscale == 0
        warning('units unspecified');
        if mean(std(D(P.Ic,ltb,ltr)))>1e-2
            guess_ind=[1 2 3];
        else
            guess_ind=[3 2 1];
        end
        msg_str  = sprintf('Units of EEG are %s ? (rms=%3.2e)',allunits(guess_ind(1),:),mean(std(D(P.Ic,ltb,ltr))));
        dip_ch   = sprintf('%s|%s|%s',allunits(guess_ind(1),:),allunits(guess_ind(2),:),allunits(guess_ind(3),:));
        dip_val  = [1,2,3];
        def_opt  = 1;
        unitind  = spm_input(msg_str,2,'b',dip_ch,dip_val,def_opt);
        allunits(guess_ind(unitind),:)
        D        = units(D, 1:Neegchans, allunits(guess_ind(unitind),:));
        EEGscale = allscales(guess_ind(unitind));
        D.save; % Save the new units
    end % if EEGscale==0
    
end % if eeg data

% data and principal modes
%--------------------------------------------------------------------------
Nmod    = 2;                              % maximum nnumber of modes
ntr     = numel(ltr);
ntb     = numel(ltb);
dat_y   = D(P.Ic,ltb,ltr)*EEGscale;
y       = reshape(dat_y,Neegchans,ntb*ntr);
y       = spm_detrend(y);
y       = y/std(y(:));

[u,s,v] = spm_svd(y,exp(-4));
i       = 1:min(size(u,2),Nmod);
y       = u(:,i)*s(i,i);
Nmod    = size(y,2);


% Other bits of the P structure, apart for priors and #dipoles
%==========================================================================
P.ltr = ltr;
P.Nc  = length(P.Ic);

% Deal with dipoles number and priors
%==========================================================================
Ndip  = 16;
vert  = D.inv{D.val}.mesh.tess_mni.vert;
ic    = fix(linspace(1,size(vert,1),Ndip));
vert  = vert(ic,:);
pE.p  = vert;
pE.q  = zeros(Ndip,3,Nmod);
Sp    = eye(spm_length(pE.p))*8^2;
Sq    = eye(spm_length(pE.q))/64;
pC    = blkdiag(Sp,Sq);

Lscal = 1/norm(full(spm_eeg_lgainmat(D,ic)));




% Launch inversion
%==========================================================================

% Initialise inverse field
%--------------------------------------------------------------------------
dip_c   = Ndip;
inverse = struct( ...
    'F',[], ...                % free energy
    'pst',D.time, ...          % all time points in data epoch
    'tb',tb, ...               % time window/bin used
    'ltb',ltb, ...             % list of time points used
    'ltr',ltr, ...             % list of trial types used
    'n_seeds',length(ltr), ... % using this field for multiple reconstruction
    'n_dip',dip_c, ...         % number of dipoles used
    'loc',[], ...              % loc of dip (3 x n_dip)
    'j',[], ...                % dipole(s) orient/ampl, in 1 column
    'cov_loc',[], ...          % cov matrix of source location
    'cov_j',[], ...            % cov matrix of source orient/ampl
    'Mtb',1, ...               % ind of max EEG power in time series, 1 as only 1 tb.
    'exitflag',[], ...         % Converged (1) or not (0)
    'P',[]);                   % save all kaboodle too.

Y.y  = y;                  % data: spatial profile of temporal modes


% set up figures
%----------------------------------------------------------------------
P.handles.hfig                 = spm_figure('GetWin','Graphics');
spm_clf(P.handles.hfig)
P.handles.SPMdefaults.col      = get(P.handles.hfig,'colormap');
P.handles.SPMdefaults.renderer = get(P.handles.hfig,'renderer');
set(P.handles.hfig,'userdata',P)
Niter = 8;
for j = 1:Niter
    
    % get lead fields
    %------------------------------------------------------------------
    M.pE    = pE;              % prior expectations of dipole parameters
    M.pC    = pC;              % prior covariance of dipole parameters
    M.hE    = 2;               % expected log precision of data
    M.hC    = 1/128;           % variability of the above precision
    M.Lscal = Lscal;           % lead field scaling
    M.Setup = P;               % pass volume conductor and sensor locations on
    M.IS    = 'spm_eeg_wrap_dipfit_vbecd_new';
    M.Nmax  = 2;
    
    % invert
    %----------------------------------------------------------------------
    [Ep,Cp,Eh,F] = spm_nlsi_GN(M,[],Y);
    
    
    % Bayesian model reduction to remove redundant dipoles
    %----------------------------------------------------------------------
    GF = [];
    while Ndip > 4
        G     = zeros(Ndip,1);
        for i = 1:Ndip
            rE          = spm_zeros(pE);
            rE.q(i,:,:) = 1;
            R           = diag(spm_vec(rE) == 0);
            rE          = R*spm_vec(pE);
            rC          = R*pC*R;
            G(i)        = spm_log_evidence(Ep,Cp,pE,pC,rE,rC);
        end
        
        % remove redundant dipole
        %------------------------------------------------------------------
        [q,i]       = max(G);
        rE          = spm_zeros(pE);
        rE.q(i,:,:) = 1;
        R           = diag(spm_vec(rE) == 0);
        rE          = R*spm_vec(pE);
        rC          = R*pC*R;
        
        % updated priors and posteriors
        %------------------------------------------------------------------
        [G,Ep,Cp]   = spm_log_evidence(Ep,Cp,pE,pC,rE,rC);
        pE          = spm_unvec(rE,Ep);
        pC          = rC;
        
        % remove redundant dipole
        %------------------------------------------------------------------
        rE          = spm_zeros(pE);
        rE.p(i,:)   = 1;
        rE.q(i,:,:) = 1;
        r           = find(spm_vec(rE) == 0);
        pE.p(i,:)   = [];
        pE.q(i,:,:) = [];
        Ep.p(i,:)   = [];
        Ep.q(i,:,:) = [];
        pC          = pC(r,r);
        Cp          = Cp(r,r);
        
        Ndip = Ndip - 1;
        GF(end + 1) = G;
        
    end
    
    ypost       = spm_eeg_wrap_dipfit_vbecd_new(Ep,M);
    Nparam      = spm_length(Ep.p);
    P.F         = F;
    P.y         = Y.y(:,1);
    P.ypost     = ypost(:,1);
    P.post_mu_s = Ep.p';
    P.post_mu_w = Ep.q(:,:,1)';
    P.post_S_s  = Cp(1:Nparam,1:Nparam);
    P.post_S_w  = Cp(Nparam+1:end,Nparam+1:end);
    
    varresids   = var(P.y - P.ypost);
    pov         = 100*(1-varresids(j)/var(P.y)); % percent variance explained

    
    % display
    %--------------------------------------------------------------
    megloc    = P.post_mu_s; % loc of dip (3 x n_dip)
    mniloc    = D.inv{val}.datareg(modind).toMNI*[megloc;ones(1,size(megloc,2))]; % actual MNI location (with scaling)
    megmom    = P.post_mu_w; % moments of dip (3 x n_dip)
    megposvar = spm_unvec(diag(P.post_S_s),P.post_mu_s); % estimate of positional uncertainty in three principal axes
    mnimom    = orM1*megmom; % convert moments into mni coordinates through a rotation (no scaling or translation)
    mniposvar = (orM1*sqrt(megposvar)).^2; % convert pos variance into approx mni space by switching axes
    
    displayVBupdate2(P.y,mnimom,mniloc(1:3,:),mniposvar,P,P.ypost);
    
end; % for  j

% Get the results out
%----------------------------------------------------------------------
inverse.pst      = tb*1e3;
inverse.F        = P.F; % free energy
inverse.mniloc   = mniloc(1:3,:);
inverse.loc      = megloc;

inverse.j        = P.post_mu_w; % dipole(s) orient/ampl, in 1 column in meg space
inverse.jmni     = reshape(mni_w,1,numel(mni_w))'; % dipole(s) orient/ampl in mni space
inverse.cov_loc  = P.post_S_s; % cov matrix of source location
inverse.cov_j    = P.post_S_w; % cov matrix of source orient/ampl
inverse.exitflag = 1; % Converged (1) or not (0)
inverse.P        = P; % save all kaboodle too.


% Save results and display
%--------------------------------------------------------------------------
D.inv{val}.inverse = inverse;
save(D)

return


function [P] = displayVBupdate2(y,mu_w,mu_s,diagS_s,P,yHat)
% yHat is estimate of y based on dipole position
%--------------------------------------------------------------------------

% plot dipoles
%----------------------------------------------------------------------
try
    opt.ParentAxes = P.handles.axesECD;
    opt.hfig       = P.handles.hfig;
    opt.handles.hp = P.handles.hp;
    opt.handles.hq = P.handles.hq;
    opt.handles.hs = P.handles.hs;
    opt.handles.ht = P.handles.ht;
    opt.query = 'replace';
catch
    P.handles.axesECD = axes(...
        'parent',P.handles.hfig,...
        'Position',[0.13 0.55 0.775 0.4],...
        'hittest','off',...
        'visible','off');
    opt.ParentAxes = P.handles.axesECD;
    opt.hfig = P.handles.hfig;
end
w = reshape(mu_w,3,[]);
s = reshape(mu_s,3,[]);
[out] = spm_eeg_displayECD(s,w,diagS_s,[],opt);

P.handles.hp = out.handles.hp;
P.handles.hq = out.handles.hq;
P.handles.hs = out.handles.hs;
P.handles.ht = out.handles.ht;


% plot data and predicted data
%--------------------------------------------------------------------------
pos          = P.forward.sens.prj;
ChanLabel    = P.channels;
in.f         = P.handles.hfig;
in.noButtons = 1;
try
    P.handles.axesY;
catch
    figure(P.handles.hfig)
    P.handles.axesY = axes(...
        'Position',[0.02 0.3 0.3 0.2],...
        'hittest','off');
    in.ParentAxes = P.handles.axesY;
    spm_eeg_plotScalpData(y,pos,ChanLabel,in);
    title(P.handles.axesY,'measured data')
end


miY = min([yHat;y]);
maY = max([yHat;y]);
try
    P.handles.axesYhat;
    d    = get(P.handles.axesYhat,'userdata');
    yHat = yHat(d.goodChannels);
    clim = [min(yHat(:))-( max(yHat(:))-min(yHat(:)) )/63,max(yHat(:))];
    ZI   = griddata(...
        d.interp.pos(1,:),d.interp.pos(2,:),full(double(yHat)),...
        d.interp.XI,d.interp.YI);
    set(d.hi,'Cdata',flipud(ZI));
    caxis(P.handles.axesYhat,clim);
    delete(d.hc)
    [~,d.hc] = contour(P.handles.axesYhat,flipud(ZI),...
        'linecolor',0.5.*ones(3,1));
    set(P.handles.axesYhat,...
        'userdata',d);
catch
    figure(P.handles.hfig)
    P.handles.axesYhat = axes(...
        'Position',[0.37 0.3 0.3 0.2],...
        'hittest','off');
    in.ParentAxes = P.handles.axesYhat;
    spm_eeg_plotScalpData(yHat,pos,ChanLabel,in);
    title(P.handles.axesYhat,'predicted data')
end
try
    P.handles.axesYhatY;
catch
    figure(P.handles.hfig)
    P.handles.axesYhatY = axes(...
        'Position',[0.72 0.3 0.25 0.2],...
        'NextPlot','replace',...
        'box','on');
end
plot(P.handles.axesYhatY,y,yHat,'.')
set(P.handles.axesYhatY,...
    'nextplot','add')
plot(P.handles.axesYhatY,[miY;maY],[miY;maY],'r')
set(P.handles.axesYhatY,...
    'nextplot','replace')
title(P.handles.axesYhatY,'predicted vs measured data')
axis(P.handles.axesYhatY,'square','tight')
grid(P.handles.axesYhatY,'on')






% plot precision hyperparameters
%----------------------------------------------------------------------
try
    P.handles.axesVar1;
catch
    figure(P.handles.hfig)
    P.handles.axesVar1 = axes(...
        'Position',[0.05 0.05 0.25 0.2],...
        'NextPlot','replace',...
        'box','on');
end

title(P.handles.axesVar1,'Free energy ')
axis(P.handles.axesVar1,'square');
set(P.handles.axesVar1,'Ylimmode','auto'); %,'tight')


try
    P.handles.axesVar2;
catch
    figure(P.handles.hfig)
    P.handles.axesVar2 = axes(...
        'Position',[0.37 0.05 0.25 0.2],...
        'NextPlot','replace',...
        'box','on');
end

title(P.handles.axesVar2,'Percent variance explained');
axis(P.handles.axesVar2,'square');


try
    P.handles.axesVar3;
catch
    figure(P.handles.hfig)
    P.handles.axesVar3 = axes(...
        'Position',[0.72 0.05 0.25 0.2],...
        'NextPlot','replace',...
        'box','on');
end







