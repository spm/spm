function D = spm_eeg_inv_vb_ecd_gui(D,val)
% GUI function for variational Bayesian ECD inversion
%
% Fills in the following fields of the inverse structure:
% inverse = struct( ...
%     'F',            % free energies as dipoles are removed
%     'pst',          % all time points in data epoch
%     'tb',           % time window/bin used
%     'ltb',          % list of time points used
%     'ltr',          % list of trial types used
%     'n_dip',        % number of dipoles used
%     'Lecd',         % dipole lead fields
%     'loc',          % loc of dipoles (n_dip x 3)
%     'exitflag',     % Converged (1) or not (0)
%     'P'             % forward model
%
% In brief, this routine:
% - load the necessary data, if not provided,
% - fill in all the necessary bits for the VB-ECD inversion routine,
% - launch variational Bayesian model inversion,
% - eliminates redundant dipoles using Bayesian model reduction,
% - displays the results.
%
% This routine provides a Bayes optimal solution to the ECD problem. It
% finesses the nonlinear inversion problem by starting with a large number
% of dipoles (on the cortical surface). It then fits the principal spatial
% modes of the data over a specified peristimulus time window using fixed
% dipole orientations. Finally, it uses Bayesian model reduction to
% eliminate the least likely dipoles, until the specified number of dipoles
% is obtained.
%
% The purpose of this routine is to find the location of a small number of
% dipoles that accurately explain fluctuations in activity over
% peristimulus time. It is anticipated that the moments of the dipoles will
% be estimated as needed using a standard pseudo-inverse (ordinary least
% squares) estimator - should it be required. examples of this are provided
% during the presentation of the results below.
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_eeg_inv_vb_ecd_gui.m 6966 2016-12-09 10:26:26Z guillaume $


% Load data, if necessary
%==========================================================================
if nargin < 1
    D = spm_eeg_load;
end

% Check if the forward model was prepared & handle the other info bits
%==========================================================================
if ~isfield(D,'inv')
    error('Data must have been prepared for inversion procedure.');
end

% check index provided
%--------------------------------------------------------------------------
if nargin == 2
    if val > length(D.inv)
        val   = length(D.inv);
        D.val = val;
    end
else
    if isfield(D,'val')
        val   = D.val;
    else
        val   = length(D.inv);
        D.val = val;
    end
end

% Use val to define which is the "current" inv{} to use
% If no inverse solution already calculated (field 'inverse' doesn't exist)
% use that inv{}. Otherwise create a new one from the previous
%--------------------------------------------------------------------------
if isfield(D.inv{val},'inverse')
    Ninv          = length(D.inv);
    D.inv{Ninv+1} = D.inv{val};
    if isfield(D.inv{Ninv + 1},'contrast')
        D.inv{Ninv + 1} = rmfield(D.inv{Ninv + 1},'contrast');
    end
    val   = Ninv + 1;
    D.val = val;
end

% Set time, date, comments & modality
%--------------------------------------------------------------------------
if ~isfield(D.inv{val},'date')
    clck = fix(clock);
    if clck(5) < 10
        clck = [num2str(clck(4)) ':0' num2str(clck(5))];
    else
        clck = [num2str(clck(4)) ':' num2str(clck(5))];
    end
    D.inv{val}.date = char(date,clck);
end

D.inv{val}.method = 'vbecd';

% Struct that collects the inputs for vbecd code
%--------------------------------------------------------------------------
P                 = [];
P.modality        = spm_eeg_modality_ui(D, 1, 1);
data              = spm_eeg_inv_get_vol_sens(D,val,[],'inv',P.modality);
P.forward.vol     = data.(P.modality(1:3)).vol;
if ischar(P.forward.vol)
    P.forward.vol = ft_read_vol(P.forward.vol);
end
P.forward.sens    = data.(P.modality(1:3)).sens;
P.forward.siunits = data.siunits;

% check channels
%--------------------------------------------------------------------------
P.Ic    = indchantype(D, P.modality, 'GOOD');
if isempty(P.Ic)
    error(['The specified modality (' P.modality ') is missing from file ' D.fname]);
else
    P.channels = D.chanlabels(P.Ic);
end

P.forward.chanunits = D.units(P.Ic);
[P.forward.vol, P.forward.sens] = ft_prepare_vol_sens( ...
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
        if all(tb >= floor(min(D.time([], 'ms')))) && all(tb<=ceil(max(D.time([], 'ms')))) && tb(1)<=tb(2)
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

% number of dipoles
%--------------------------------------------------------------------------
Sdip = spm_input('number of dipoles',2,'b',{'1','2','4','8','16','32'});
Sdip = str2double(Sdip);

% restrict source locations (see table at the end of this script)
%--------------------------------------------------------------------------
Sres = spm_input('restrict sources',3,'b',{'none','auditory','visual'});
switch Sres
    case 'auditory'
        res(1).centre = [ 56 4 10];
        res(1).radius = 36;
        res(2).centre = [-56 4 10];
        res(2).radius = 36;
        
    case 'visual'
        res(1).centre = [ 28 96 6];
        res(1).radius = 36;
        res(2).centre = [-28 96 6];
        res(2).radius = 36;
        
    otherwise
        res = [];
end


% trial type
%--------------------------------------------------------------------------
if D.ntrials > 1
    msg_tr = ['Trial type number [1 ',num2str(D.ntrials),']'];
    ltr    = spm_input(msg_tr,4,'i',num2str(1:D.ntrials));
else
    ltr    = 1;
end

% data, averaged over time window considered
%--------------------------------------------------------------------------
EEGscale = 1;

% SORT OUT EEG UNITS AND CONVERT VALUES TO VOLTS
%--------------------------------------------------------------------------
if strcmpi(P.modality,'EEG')
    allunits  = strvcat('uV','mV','V');
    allscales = [1, 1e3, 1e6];
    EEGscale  = 0;
    eegunits  = unique(D.units(D.indchantype('EEG')));
    Nchan     = numel( D.units(D.indchantype('EEG')));
    for j=1:length(allunits)
        if strcmp(deblank(allunits(j,:)),deblank(eegunits))
            EEGscale=allscales(j);
        end
    end
    
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
        unitind  = spm_input(msg_str,4,'b',dip_ch,dip_val,def_opt);
        allunits(guess_ind(unitind),:)
        D        = units(D, 1:Nchan, allunits(guess_ind(unitind),:));
        EEGscale = allscales(guess_ind(unitind));
        D.save; % Save the new units
    end
    
end % if eeg data

% data and principal modes
%--------------------------------------------------------------------------
Nmod  = 4;                              % maximum number of modes
ntr   = numel(ltr);
ntb   = numel(ltb);
dat_y = D(P.Ic,ltb,ltr)*EEGscale;
y     = reshape(dat_y,Nchan,ntb*ntr);
y     = spm_detrend(y);
y     = y/std(y(:));

[u,s] = spm_svd(y,exp(-4));
i     = 1:min(size(u,2),Nmod);
Y.y   = u(:,i)*s(i,i);
Nmod  = size(Y.y,2);


% Other bits of the P structure, apart for priors and #dipoles
%==========================================================================
P.ltr = ltr;
P.Nc  = length(P.Ic);

% get dipole locations and lead fields
%==========================================================================
Ndip  = 64;
vert  = D.inv{D.val}.mesh.tess_mni.vert;

% restrict sources
%--------------------------------------------------------------------------
r     = zeros(size(vert,1),1);
for i = 1:length(res)
    j = vert - ones(size(vert,1),1)*res(i).centre;
    r = r | sum(j.^2,2) < res(i).radius^2;
end
vert  = vert(r,:);
is    = fix(linspace(1,size(vert,1),Ndip));
pos   = vert(is,:);

% forward model
%--------------------------------------------------------------------------
L     = spm_eeg_lgainmat(D,is);
L     = full(L/std(L(:)));

% Launch inversion
%==========================================================================

% priors
%--------------------------------------------------------------------------
pE     = zeros(Ndip,Nmod);

% symmetries rooms constraints
%--------------------------------------------------------------------------
for i = 1:Ndip
    for j = 1:Ndip
        d(i,j)  = sum((pos(i,:) - pos(j,:)).^2);
        h(i,j)  = sum([pos(i,1) + pos(j,1), pos(i,2) - pos(j,2), pos(i,3) - pos(j,3)].^2);
        pC(i,j) = exp(-min(d(i,j),h(i,j))/128);
    end
end
pC     = kron(eye(Nmod),pC);


% Set up model
%--------------------------------------------------------------------------
M.pE   = pE;              % prior expectations of dipole parameters
M.pC   = pC;              % prior covariance of dipole parameters
M.hE   = 2;               % expected log precision of data
M.hC   = 1/128;           % variability of the above precision
M.L    = L;               % lead field scaling
M.IS   = @(P,M,U)(M.L*P);
M.Nmax = 16;

% invert
%--------------------------------------------------------------------------
[Ep,Cp] = spm_nlsi_GN(M,[],Y);


% Bayesian model reduction to remove redundant dipoles
%--------------------------------------------------------------------------
F = [];
N = [];
while Ndip > Sdip
    G     = zeros(Ndip,1);
    for r = 1:Ndip
        rE      = spm_zeros(pE);
        rE(r,:) = 1;
        R       = diag(~spm_vec(rE));
        rE      = R*spm_vec(pE);
        rC      = R*pC*R';
        G(r)    = spm_log_evidence(Ep,Cp,pE,pC,rE,rC);
    end
    
    % select the priors of the most likely model
    %----------------------------------------------------------------------
    [~,i] = sort(-G);
    if Ndip > 128
        r = G < G(i(64));
    elseif Ndip > 32
        r = G < G(i(8));
    else
        r = G < G(i(1));
    end
    rE        = spm_zeros(pE);
    rE(r,:)   = 1;
    R         = diag(spm_vec(rE));
    rE        = R*spm_vec(pE);
    rC        = R*pC*R';
    
    % update priors and posteriors
    %----------------------------------------------------------------------
    [G,Ep,Cp] = spm_log_evidence_reduce(Ep,Cp,pE,pC,rE,rC);
    pE        = spm_unvec(rE,Ep);
    pC        = rC;
    
    % remove redundant dipole from reduced priors and posterior is
    %----------------------------------------------------------------------
    r         = find(r);
    R         = R(any(R),:);
    pE        = pE(r,:);
    Ep        = Ep(r,:);
    L         = L(:,r);
    pos       = pos(r,:);
    pC        = R*pC*R';
    Cp        = R*Cp*R';
    
    % record changes in free energy
    %----------------------------------------------------------------------
    Ndip       = numel(r);
    F(end + 1) = G;
    N(end + 1) = Ndip;
    fprintf('number of sources %i (dF = %0.2f)\n',Ndip,G)
    
end

% accumulate changes in free energy
%--------------------------------------------------------------------------
F = cumsum(F);


% set up figure and show results
%==========================================================================
P.handles.hfig                 = spm_figure('GetWin','Graphics');
spm_clf(P.handles.hfig)
P.handles.SPMdefaults.col      = get(P.handles.hfig,'colormap');
P.handles.SPMdefaults.renderer = get(P.handles.hfig,'renderer');
set(P.handles.hfig,'userdata',P)

% show dipoles
%--------------------------------------------------------------------------
opt.ParentAxes = subplot(2,2,1);
opt.hfig = P.handles.hfig;
spm_eeg_displayECD(pos',Ep(:,1)',8,[],opt);


% show free energy following Bayesian model reduction
%--------------------------------------------------------------------------
subplot(4,2,2), plot(N,F), axis square, spm_axis tight
xlabel('number of dipoles remaining'), 
ylabel('change in log evidence'), title('Bayesian model reduction')


% on source activity
%--------------------------------------------------------------------------
pst   = D.time(ltb)*1000;
for i = 1:ntr
    J = pinv(L)*dat_y(:,:,i);
    subplot(4,2,4), plot(pst,J), hold on
end
xlabel('number of dipoles remaining'), axis square, spm_axis tight, hold off
ylabel('change in log evidence'), title('Bayesian model reduction')



% plot observed and predicted primary mode
%--------------------------------------------------------------------------
Yobs  = Y.y(:,1);
Ypred = L*(pinv(L)*Yobs);

in.f          = P.handles.hfig;
in.noButtons  = 1;
in.ParentAxes = subplot(4,2,5);
spm_eeg_plotScalpData(Yobs,P.forward.sens.prj,P.channels,in);
title('measured mode')

in.ParentAxes = subplot(4,2,6);
spm_eeg_plotScalpData(Ypred,P.forward.sens.prj,P.channels,in);
title('predicted mode')

% plot observed and predicted responses
%--------------------------------------------------------------------------
Yobs  = y;
Ypred = L*(pinv(L)*Yobs);
pov   = (var(Yobs(:)) - var(Yobs(:) - Ypred(:)))/var(Yobs(:));
pov   = sprintf('variance explained %0.4g',pov*100);

subplot(4,2,7); imagesc(repmat(pst,1,Nmod),P.Ic,Yobs)
xlabel('peristimulus time (ms)'), ylabel('channels'), title('measured data' )
subplot(4,2,8); imagesc(repmat(pst,1,Nmod),P.Ic,Ypred)
xlabel(pov), ylabel('channels'), title('predicted data')
    
    
% inverse field
%--------------------------------------------------------------------------
inverse = struct( ...
    'F',F, ...                 % free energies over dipoles
    'pst',D.time, ...          % all time points in data epoch
    'tb',pst, ...              % time window/bin used
    'ltb',ltb, ...             % list of time points used
    'ltr',ltr, ...             % list of trial types used
    'n_dip',Ndip, ...          % number of dipoles used
    'Lecd',L, ...              % dipole lead fields
    'loc',pos, ...             % loc of dip (n_dip x 3)
    'exitflag',1, ...          % Converged (1) or not (0)
    'n_seeds',1, ...           % one 
    'Mtb',1, ...               % one 
    'P',[]);                   % save all kaboodle too.

% Save results and display
%--------------------------------------------------------------------------
inverse.mniloc = {pos'};
D.inv{val}.inverse = inverse;
save(D)

return


% MNI coordinates of 64 regions of interest
%-------------------------------------------------------------------------
% Region Component  b X Y Z 
% Gray matter anterior intraparietal sulcus hIP1 right 1 32 60 50
% Gray matter anterior intraparietal sulcus hIP1 left 1 24 60 46
% Gray matter primary auditory cortex TE1.2 right 1 56 4 10
% Gray matter primary auditory cortex TE 1.2 left 1 56 4 10
% Gray matter visual cortex V1 BA17 right 2 28 96 6
% Gray matter visual cortex V1 BA17 left 2 28 96 6
% Gray matter Broca area BA44 right 3 44 8 54
% Gray matter Broca area BA44 left 3 38 16 50
% Gray matter primary motor cortex BA4a right 3 52 4 54 
% Gray matter primary motor cortex BA4a left 3 48 8 50
% Gray matter premotor cortex BA6 right 3 28 20 58 
% Gray matter premotor cortex BA6 left 3 24 20 58 
% Gray matter visual cortex V5 right 4 44 84 10 
% Gray matter visual cortex V5 left 4 36 86 8 
% Gray matter anterior intraparietal sulcus hIP1 right 6 52 40 48 
% Gray matter anterior intraparietal sulcus hIP1 left 6 48 38 42 
% Gray matter anterior intraparietal sulcus hIP2 right 6 66 28 30 
% Gray matter anterior intraparietal sulcus hIP2 left 6 60 32 26 
% Gray matter premotor right 6 60 8 22 
% Gray matter premotor left 6 52 4 22 
% Gray matter primary motor cortex BA4p right 6 62 14 30 
% Gray matter primary motor cortex BA4p left 6 56 20 30 
% Gray matter visual cortex V5 right 6 52 68 10 
% Gray matter visual cortex V5 left 6 48 68 10 
% Gray matter anterior intraparietal sulcus hIP1 right 8 50 60 26 
% Gray matter anterior intraparietal sulcus hIP1 left 7 48 58 30 
% Gray matter Broca area BA45 right 8 52 22 10 
% Gray matter Broca area BA45 left 7 48 24 10 
% Gray matter amygdala-centromedial group right 41 24 8 18 
% Gray matter amygdala-centromedial group left 41 20 2 20 
% Gray matter anterior intraparietal sulcus hIP1 right 14 26 44 44 
% Gray matter anterior intraparietal sulcus hIP1 left 14 28 42 40 
% Gray matter primary motor cortex BA4a right 14 10 44 46 
% Gray matter primary motor cortex BA4a left 14 8 40 46 
% Gray matter amygdala-laterobasal group right 13 38 2 38 
% Gray matter amygdala-laterobasal group left 13 38 2 34 
% Gray matter primary auditory cortex TE1.0 right 18 48 28 10 
% Gray matter primary auditory cortex TE1.0 left 18 40 28 6 
% Gray matter primary motor cortex BA4a right 19 16 24 62 
% Gray matter primary motor cortex BA4p left 19 16 28 62 
% Gray matter anterior intraparietal sulcus hIP1 right 23 52 60 44 
% Gray matter anterior intraparietal sulcus hIP1 left 23 46 60 46 
% Gray matter primary motor cortex BA4a right 24 12 16 74 
% Gray matter primary motor cortex BA4a left 24 8 14 74 
% Gray matter primary motor cortex BA4p right 24 32 28 70 
% Gray matter primary motor cortex BA4p left 24 28 28 70 
% Gray matter primary somatosensory cortex BA3a right 24 8 34 70 
% Gray matter primary somatosensory cortex BA3a left 24 8 38 70 
% Gray matter primary auditory cortex TE1.0 right 25 64 24 6 
% Gray matter primary auditory cortex TE1.0 left 25 60 28 6 
% Gray matter Broca area BA45 right 31 48 38 16 
% Gray matter Broca area BA44 left 31 40 38 22 
% Gray matter hippocampus cornu ammonis right 32 32 40 16 
% Gray matter hippocampus cornu ammonis left 32 28 44 14 
% Gray matter anterior intraparietal sulcus hIP2 right 38 20 58 22 
% Gray matter anterior intraparietal sulcus hIP1 left 38 16 60 22 
% Gray matter primary somatosensory cortex BA3a right 40 12 52 68 
% Gray matter primary somatosensory cortex BA3a left 40 8 54 62 
% Gray matter visual cortex V5 right 66 50 64 8 
% Gray matter visual cortex V5 left 66 50 64 8 
% Frontal Pole 54 0 64 18
% Superior frontal gyrus 54 0 40 50
% Frontal medial cortex 54 4 55 10
% Gray matter visual cortex V1 BA17 63 2 80 34
% a Coordinates are adapted from Kiviniemi et al.32
% b The independent component to which the coordinate is localized.
% c The region that is the interhemispheric homologue to a given region.
% AJNR Am J Neuroradiol   www.ajnr.org 1
% On-line Fig 1. Graphic illustration of 64 ROIs from On-line Table 1 used for connectivity analysis.
% 2  AJNR   www.ajnr.org
