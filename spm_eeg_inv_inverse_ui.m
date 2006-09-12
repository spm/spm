function D = spm_eeg_inv_inverse_ui(D)

%=======================================================================
% Inverse Solution user-interface routine
% commands the ReML inverse computation for either EEG or MEG data
% to reconstruct the sources of either evoked or induced activity.
%
% FORMAT D = spm_eeg_inv_inverse_ui(S)
% Input:
% S		    - input data struct (optional)
% Output:
% D			- same data struct including the inverse solution files and variables
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_inverse_ui.m 621 2006-09-12 17:22:42Z karl $

spm_defaults

try
    D;
catch
    D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
    D = spm_eeg_ldata(D);
end

try
    val = D.val;
catch
    val = length(D.inv);
end


% Type of analysis
if isempty(D.inv{val}.inverse.activity)
    activity = spm_input('Type of analysis','+1','evoked|induced|both',[1 2 3]);
    switch activity
        case 1 % estimate the evoked activity (require averaged data)
            if D.events.Ntypes ~= D.Nevents
                error(sprintf('Load averaged data instead\n'));
            else
                D.inv{val}.inverse.activity = 'evoked';
            end
        case 2 % estimate the induced activity (require single trial data)
            if length(D.events.code) ~= D.Nevents
                error(sprintf('Wrong data structure\n'));
            else
                D.inv{val}.inverse.activity = 'induced';
            end
        case 3 % estimate both the evoked and induced activity (requires single trial data)
            if length(D.events.code) ~= D.Nevents
                error(sprintf('Wrong data structure\n'));
            else
                D.inv{val}.inverse.activity = 'evoked & induced';
            end
    end
end

% Contrast of interest between conditions
if isempty(D.inv{val}.inverse.contrast)
    if D.events.Ntypes > 1
        contrast = spm_input(sprintf('1 x %i Contrast vector',D.events.Ntypes),'+1');
        while (length(contrast) ~= D.events.Ntypes) | ~isnumeric(contrast)
            warndlg({'Wrong contrast vector format','please try again'});
            contrast = spm_input(sprintf('1 x %i Contrast vector',D.events.Ntypes),'0');
        end
        if strcmp(D.inv{val}.inverse.activity,'induced') & length(setdiff(contrast,zeros(size(contrast)))) > 1
            warndlg({'Select 1 trial type only for induced responses','Please try again'});
            return
        end
    else
        disp('There is only one condition');
        contrast = [1];
    end
    D.inv{val}.inverse.contrast = contrast;
    clear contrast
end

% Time window of interest
if isempty(D.inv{val}.inverse.woi)
    Tstart = -D.events.start/(D.Radc/1000);
    Tstop  =  D.events.stop /(D.Radc/1000);
    woi = spm_input('Time window of interest (ms)','+1','r',round([Tstart Tstop]));
    if woi(1) > woi(2)
        error(sprintf('Wrong entry!\n'));
    end
    Sstart = round(woi(1)*(D.Radc/1000));
    if (woi(1) < 0) & (abs(Sstart) > D.events.start)
        Sstart = -D.events.start;
    end
    Sstop = round(woi(2)*(D.Radc/1000));
    if Sstop > D.events.stop
        Sstop = D.events.stop;
    end
    D.inv{val}.inverse.woi = round([Sstart Sstop]/(D.Radc/1000));
    clear woi Tstart Tstop Sstart Sstop
end


if strcmp(D.inv{val}.inverse.activity,'induced') | strcmp(D.inv{val}.inverse.activity,'evoked & induced')
    % Frequency band of interest
    if isempty(D.inv{val}.inverse.fboi)
        Flow  = 8;
        Fhigh = 12;
        fboi  = spm_input('Frequency band of interest (Hz)','+1','r',[Flow Fhigh]);
        if fboi(1) <= 0 | fboi(2) <= 0 | fboi(1) > fboi(2)
            error(sprintf('Wrong entry!\n'));
        end
        D.inv{val}.inverse.fboi = [Flow Fhigh];
    end
end


% Source space reduction
% Can be specified as a ratio (between 0 and 1) wrt the initial number of
% nodes / Or as the integer value (greater than 1) specifying the
% new number of nodes (less than the initial number)

% switch this off for general users
%--------------------------------------------------------------------------
D.inv{val}.inverse.dim = D.inv{val}.mesh.Ctx_Nv;
if isempty(D.inv{val}.inverse.dim)
    Rfactor = spm_input('Source space reduction','+1','r',1);
    if Rfactor > D.inv{val}.mesh.Ctx_Nv | Rfactor <= 0
        disp('Impossible value');
        return
    elseif Rfactor <= 1
        Nnodes = ceil(Rfactor*D.inv{val}.mesh.Ctx_Nv);
    else
        Nnodes = Rfactor;
    end
    D.inv{val}.inverse.dim = Nnodes;
    clear Rfactor Nnodes
end

% Define and compute covariance components
% Sensor space (level 1)
% evoked
%--------------------------------------------------------------------------
if strcmp(D.inv{val}.inverse.activity,'evoked')
    str = {'i.i.d (1st level)',...
           '& baseline noise',...
           '& anti-averging'};
    Vs  = spm_input('sensor-level components','+1','m',str);
    if     Vs == 1
        Vsens = [1 0 0];
    elseif Vs == 2
        Vsens = [1 0 1];
    elseif Vs == 3
        Vsens = [1 1 1];
    else
        Vsens = [1 0 0];
    end
    
% evoked induced
%--------------------------------------------------------------------------
else
    Vsens = [1 0 0];
end

% Source space (level 2)
% evoked
%--------------------------------------------------------------------------
if strcmp(D.inv{val}.inverse.activity,'evoked')
    str = {'Smoothness (2nd level)',...
           '& Minimum Norm'};
    Vs  = spm_input('source-level components','+1','m',str);
    if     Vs == 1
        Vsour = [1 0 0];
    elseif Vs == 2
        Vsour = [1 1 0];
    else
        Vsour = [1 0 0];
    end
    
    if spm_input('& user-specifed source prior','+1','yes|no',[1 0]);
        Vsour = [Vsour 1];
    end
    
% evoked
%--------------------------------------------------------------------------
else
    str = {'Smoothness (2nd level)',...
            '& Minimum Norm'};
    Vs  = spm_input('source-level components','+1','m',str);
    if     Vs == 1
        Vsour = [1 0 0];
    elseif Vs == 2
        Vsour = [1 1 0];
    else
        Vsour = [1 0 0];
    end
    
    if spm_input('& user-specified source prior','+1','yes|no',[1 0]);
        Vsour = [Vsour 1];
    end
end

if length(Vsour) > 3
    P = spm_select(1, '.mat', 'Select prior file');
    D.inv{val}.inverse.priors.level2{4}.filename = P;
    l2_covar = load(P);
    name     = fieldnames(l2_covar);
    D.inv{val}.inverse.priors.level2{4}.label  = name{1};
    D.inv{val}.inverse.priors.level2{4}.status = 'yes';
    clear l2_covar name
end

% Proceed with the empirical Bayesian inference
D = spm_eeg_inv_inverse(D,Vsens,Vsour);

save(fullfile(D.path,D.fname),'D');

return