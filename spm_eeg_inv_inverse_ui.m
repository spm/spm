function D = spm_eeg_inv_inverse_ui(S)

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
% $Id: spm_eeg_inv_inverse_ui.m 338 2005-11-30 13:55:04Z guillaume $

spm_defaults


try
    D = S; 
    clear S
catch
    D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
	D = spm_eeg_ldata(D);
end


val = length(D.inv);


if val > 1
    Atype = spm_input('Set new variance comp. only','+1','yes|no',[1 0]);
else
    Atype = 0;
end


if Atype == 1

    % Copy fields from previous analysis
    D = spm_eeg_inv_copyfields(D,[0 0 0 1],0);
    
else

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
            disp(['There are ' num2str(D.events.Ntypes) ' event types']);
            contrast = spm_input('Contrast vector','+1');
            if (length(contrast) ~= D.events.Ntypes) | ~isnumeric(contrast)
                error(sprintf('Wrong contrast vector format\n'));
            end
            if strcmp(D.inv{val}.inverse.activity,'induced') & length(setdiff(contrast,zeros(size(contrast)))) > 1
                error(sprintf('Select 1 trial type only\n'));
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
        Tstop  = D.events.stop/(D.Radc/1000);
        woi = spm_input('Window of interest (ms)','+1','r',round([Tstart Tstop]));
        if woi(1) > woi(2)
            error(sprintf('Wrong entry!\n'));
        end
        Sstart = round(woi(1)*(D.Radc/1000));
        if (woi(1) < 0) & (abs(Sstart) > D.events.start)
            Sstart = -D.events.start;
        end
        if woi(2) <= 0
            error(sprintf('Wrong entry!\n'));
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
            Flow   = 8;
            Fhigh  = 12;
            fboi = spm_input('Frequecny band of interest (Hz)','+1','r',[Flow Fhigh]));
            if fboi(1) <= 0 | fboi(2) <= 0 | fboi(1) > fboi(2)
                error(sprintf('Wrong entry!\n'));
            end
            D.inv{val}.inverse.fboi = [Flow Fhigh];
        end
    end
    
end
    
% Source space reduction
% Can be specified as a ratio (between 0 and 1) wrt the initial number of
% nodes / Or as the integer value (greater than 1) specifying the
% new number of nodes (less than the initial number)
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
if strcmp(D.inv{val}.inverse.activity,'evoked')
    l1f1 = spm_input('i.i.d (1st level)','+1','Yes|No',[1 0]);
    l1f2 = spm_input('Anti-average (1st level)','+1','Yes|No',[1 0]);
    l1f3 = spm_input('Estimated noise (1st level)','+1','Yes|No',[1 0]);
    l1d  = spm_input('#other covar to load (1st level)','+1','i',0);
else
    l1f1 = spm_input('i.i.d','+1','Yes|No',[1 0]);
    l1f2 = 0;
    l1f3 = 0;
    l1d  = spm_input('#other covar to load (1st level)?','+1','i',0);
end
Vsens = [l1f1 l1f2 l1f3 ones(1,l1d)];
if sum(Vsens) == 0
    Vsens = [1 0 0]; % default (i.i.d component)
end
np = sum(Vsens(1:3));
for i = 1:l1d
    P = spm_select(1, '.mat', 'Select covar. file');
    D.inv{val}.inverse.priors.level1{3+i}.filename = P;
    l1_covar = load(P);
    name = fieldnames(l1_covar);
    D.inv{val}.inverse.priors.level1{3+i}.label = name{1};    
    D.inv{val}.inverse.priors.level1{3+i}.status = 'yes';
    clear l1_covar name
end
       
% Source space (level 2)
if D.inv{val}.inverse.activity == 'evoked'
    l2f1 = spm_input('Smoothness prior (2nd level)','+1','Yes|No',[1 0]);
    l2f2 = spm_input('Minimum Norm (2nd level)','+1','Yes|No',[1 0]);
    l2f3 = spm_input('Multi. Source Preloc. (2nd level)','+1','Yes|No',[1 0]);
    l2d  = spm_input('#other prior to load (2nd level)','+1','i',0);
else
    l2f1 = spm_input('Smoothness prior','+1','Yes|No',[1 0]);
    l2f2 = spm_input('Minimum Norm','+1','Yes|No',[1 0]);
    l2f3 = 0;
    l2d  = spm_input('#other prior to load (2nd level)','+1','i',0);
end
Vsour = [l2f1 l2f2 l2f3 ones(1,l2d)];
if sum(Vsour) == 0
    Vsour = [1 0 0]; % default (smoothness prior)
end
np = sum(Vsour(1:3));
for i = 1:l2d
    P = spm_select(1, '.mat', 'Select prior file');
    D.inv{val}.inverse.priors.level2{3+i}.filename = P;
    l2_covar = load(P);
    name = fieldnames(l2_covar);
    D.inv{val}.inverse.priors.level2{3+i}.label = name{1};
    D.inv{val}.inverse.priors.level2{3+i}.status = 'yes';
    clear l2_covar name
end

     
% Proceed with the empirical Bayesian inference
D = spm_eeg_inv_inverse(D,Vsens,Vsour);

save(fullfile(D.path,D.fname),'D');

return