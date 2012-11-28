function DCM = spm_dcm_phase_data(DCM)
% Get data for DCM for phase coupling
% FORMAT [DCM] = spm_dcm_phase_data(DCM)
%
% DCM    -  DCM structure
%
% Requires/requests:
%
%    DCM.xY.Dfile   - M/EEG data filename
%    DCM.Lpos       - Matrix of source locations
%    DCM.options.trials - To select particular trials (otherwise all selected)
%    DCM.options.Fdcm   - to select frequency window for analysis
%    DCM.options.Hdcm - to select time window for filtering and Hilbert transform
%                       default=[beginning of epoch, end of epoch] 
%    DCM.options.Tdcm - to select time window for phase-coupling analysis
%    DCM.options.filter_order - order of bandpass filter
%    DCM.options.Nmodes - specify sub-sampling of trials eg
%                           Nmodes=2 to use every other trial. This can
%                           speed up model fitting. Default value=1.
%
% Sets
%
%    DCM.xY.pst     - Peristimulus Time [ms] of time-frequency data
%    DCM.xY.dt      - sampling in seconds [s]
%    DCM.xY.y       - concatenated induced response over sources
%    DCM.xY.Ic      - Indices of good channels
%
%    DCM.xY.y{i}(k,l) - Phase data for i-th trial,l-th source,k-th time-bin
%
%
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_dcm_phase_data.m 5082 2012-11-28 20:25:37Z karl $

% Get data filename
%-------------------------------------------------------------------------
try
    Dfile = DCM.xY.Dfile;
catch
    errordlg('Please specify data and trials');
    error('')
end

try
    filter_order=DCM.options.filter_order;
catch
    filter_order=5;
end

% load D
%--------------------------------------------------------------------------
D = spm_eeg_load(Dfile);

pst=time(D, [], 'ms');
try
    Hdcm=DCM.options.Hdcm;
catch
    Hdcm(1)= pst(1);
    Hdcm(2)= pst(end);
end

% Check time windows are valid
Tdcm=DCM.options.Tdcm;
if Tdcm(1) < pst(1) || Tdcm(2) > pst(end)
    disp('Error in spm_dcm_phase_data: DCM.options.Tdcm is outside data time window');
    return
end
if Hdcm(1) > Tdcm(1) || Hdcm(2) < Tdcm(2)
    disp('Error in spm_dcm_phase_data: DCM.options.Hdcm is narrower than DCM.options.Tdcm')
    return
end

% DCM.options.trials (from spm_api_erp) contains list of conditions
% Get list of actual trial indices
chosen_conds=DCM.options.trials;

% Here using the otherwise unused Nmodes option to select a subset of trials
trial_step = DCM.options.Nmodes;
cond_name=condlist(D);
trials=[];
X=[];
for jj=1:length(chosen_conds)
    cname=cond_name{chosen_conds(jj)};
    new_trials= indtrial(D, cname, 'GOOD');
    
    new_trials = new_trials(1:trial_step:end);
    
    trials=[trials(:);new_trials(:)];

    X=[X;ones(length(new_trials),1)*DCM.xU.X(jj)];
end

% Change DCM.xU to accomodate these trial indices
DCM.xU.oldX=DCM.xU.X;
DCM.xU.X=X;


% indices of EEG channel (excluding bad channels)
%--------------------------------------------------------------------------
if ~isfield(DCM.xY, 'modality')
    DCM.xY.modality = spm_eeg_modality_ui(D);
end

modality = DCM.xY.modality;

if ~isfield(DCM.xY, 'Ic')
    Ic        = D.indchantype(modality,'GOOD');  
    DCM.xY.Ic = Ic;
end

Ic = DCM.xY.Ic;

% Compute indices for time window
DCM.xY.pst=time(D);
ind = D.indsample(1e-3*Hdcm);
ind = ind(1):ind(2);


% Read in trials
Ntr=length(trials);
clist=conditions(D);
for i=1:Ntr,
    k=trials(i);
    DCM.xY.y{i}=squeeze(D(Ic,ind,k))';
    DCM.xY.code{i}=clist{k};
end
DCM.xY.dt=1/fsample(D);
DCM.xY.pst=DCM.xY.pst(ind);

Nr=length(DCM.Sname);
Ntrials=length(DCM.xY.y);

if ~isequal(modality, 'LFP')
    disp(sprintf('Projecting %s data onto source locations ...',DCM.xY.modality));
    %------------------------------------------------------------------
    try
        pos = DCM.Lpos;
    catch
        pos = DCM.M.dipfit.Lpos;
    end
    Ng     = 3;
    G.L    = kron(ones(1,Nr),speye(Ng,Ng));
    G.Lpos = kron(pos,ones(1,Ng));
    L      = spm_erp_L(G,DCM.M.dipfit);
    MAP    = pinv(L);
    [Ntime,Nchannels]=size(DCM.xY.y{1});
    for n=1:Ntrials,
        % Get source time series for all orthogonal directions
        DCM.xY.y{n}=DCM.xY.y{n}*MAP';
    end

    % Get max variance orientation for each source
    for s=1:Nr,
        ind=[1:3]+(s-1)*3;
        y=[];
        for n=1:Ntrials,
            y=[y;DCM.xY.y{n}(:,ind)];
        end
        [uu,ss,vv]=svd(y,0);
        for n=1:Ntrials,
            region{n}(:,s)=DCM.xY.y{n}(:,ind)*vv(:,1);
        end
    end
else
    disp('Using data from LFP channels');
    region=DCM.xY.y;
end

DCM.xY.y = [];
% Get instantaneous phase
disp('Filter and compute instantaneous phase ...');
for n=1:Ntrials,
    for c=1:Nr,
        xr=region{n}(:,c);
        
        % Filtering
        xr = ft_preproc_bandpassfilter(xr', fsample(D), DCM.options.Fdcm, filter_order);

        hx=spm_hilbert(xr);
        DCM.xY.y{n}(:,c)=unwrap(double(angle(hx)));
    end
end
disp('Source extraction complete ...');

% Now crop data into time window of interest for phase coupling analysis
tcr=1e-3*Tdcm;
crind=find(DCM.xY.pst >= tcr(1) & DCM.xY.pst <= tcr(2));
DCM.xY.pst=DCM.xY.pst(crind);
for n=1:Ntrials,
    DCM.xY.y{n}=DCM.xY.y{n}(crind,:);
end