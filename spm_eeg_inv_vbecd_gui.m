function D = spm_eeg_inv_vbecd_gui(D,val)
% GUI function for Bayesian ECD inversion
% - load the necessary data, if not provided
% - fill in all the necessary bits for the VB-ECD inversion routine,
% - launch the B_ECD routine, aka. spm_eeg_inv_vbecd
% - displays the results.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% 
% $Id: spm_eeg_inv_vbecd_gui.m 7702 2019-11-22 11:32:26Z guillaume $

%%
% Load data, if necessary
%==========
if nargin<1
    D = spm_eeg_load;
end

DISTTHRESH=120; %% mm distance from centre of brain classified as inside the head

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
    D.inv{val}.date = strvcat(date,clck);
end

if ~isfield(D.inv{val}, 'comment'), 
   D.inv{val}.comment = {spm_input('Comment/Label for this analysis', '+1', 's')};
end

D.inv{val}.method = 'vbecd';

%% Struct that collects the inputs for vbecd code
P = [];


P.modality = spm_eeg_modality_ui(D, 1, 1);

for m = 1:numel(D.inv{val}.datareg)
    if strncmp(D.inv{val}.datareg(m).modality, P.modality, 3)
        modind = m;
    end
end

data = spm_eeg_inv_get_vol_sens(D, val, [], 'inv', P.modality);

P.forward.vol     = data.(P.modality(1:3)).vol;
if ischar(P.forward.vol)
    P.forward.vol = ft_read_headmodel(P.forward.vol);
end

P.forward.sens    = data.(P.modality(1:3)).sens;
P.forward.siunits = data.siunits;
M1 = data.transforms.toMNI;
[U, L, V] = svd(M1(1:3, 1:3));
orM1(1:3,1:3) =U*V'; %% for switching orientation between meg and mni space

P.Ic = indchantype(D, P.modality, 'GOOD');

if isempty(P.Ic)
    error(['The specified modality (' P.modality ') is missing from file ' D.fname]);
else
    P.channels = D.chanlabels(P.Ic);
end

P.forward.chanunits = D.units(P.Ic);

[P.forward.vol, P.forward.sens] =  ft_prepare_vol_sens( ...
    P.forward.vol, P.forward.sens, 'channel', P.channels);


P.forward.sens.prj = D.coor2D(P.Ic);

%% 
% Deal with data
%===============

% time bin or time window
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
if length(tb)==1
    [kk,ltb] = min(abs(D.time([], 'ms')-tb)); % round to nearest time bin
else
    [kk,ltb(1)] = min(abs(D.time([], 'ms')-tb(1)));  % round to nearest time bin
    [kk,ltb(2)] = min(abs(D.time([], 'ms')-tb(2)));
    ltb = ltb(1):ltb(2); % list of time bins 'tb' to use
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


EEGscale=1;

%% SORT OUT EEG UNITS AND CONVERT VALUES TO VOLTS
if strcmpi(P.modality,'EEG')
    allunits  = strvcat('uV','mV','V');   
    allscales = [1, 1e3, 1e6]; %% 
    EEGscale  = 0;
    eegunits  = unique(D.units(D.indchantype('EEG')));
    Neegchans = numel(D.units(D.indchantype('EEG')));
    for j=1:length(allunits)
        if strcmp(deblank(allunits(j,:)),deblank(eegunits))
            EEGscale=allscales(j);
        end % if
    end % for j
    
    if EEGscale==0
        warning('units unspecified');
        if mean(std(D(P.Ic,ltb,ltr)))>1e-2
            guess_ind=[1 2 3];
        else
            guess_ind=[3 2 1];
        end
        msg_str = sprintf('Units of EEG are %s ? (rms=%3.2e)',allunits(guess_ind(1),:),mean(std(D(P.Ic,ltb,ltr))));
        dip_ch  = sprintf('%s|%s|%s',allunits(guess_ind(1),:),allunits(guess_ind(2),:),allunits(guess_ind(3),:));
        dip_val = [1,2,3];
        def_opt = 1;
        unitind = spm_input(msg_str,2,'b',dip_ch,dip_val,def_opt);
        %ans=spm_input(msg_str,1,'s','yes');
        allunits(guess_ind(unitind),:)
        D = units(D, 1:Neegchans, allunits(guess_ind(unitind),:));
        EEGscale = allscales(guess_ind(unitind));
        D.save; %% Save the new units
    end %if EEGscale==0
    
end % if eeg data


dat_y = squeeze(mean(D(P.Ic,ltb,ltr)*EEGscale,2));


%%
% Other bits of the P structure, apart for priors and #dipoles
%==============================

P.ltr          = ltr;
P.Nc           = length(P.Ic);


% Deal with dipoles number and priors
%====================================
dip_q = 0; % number of dipole 'elements' added (single or pair)
dip_c = 0; % total number of dipoles in the model
adding_dips = 1;
clear dip_pr



priorlocvardefault=[100, 100, 100]; %% location variance default in mm
nopriorlocvardefault=[80*80, 80*80, 80*80];
nopriormomvardefault=[10, 10, 10]*100; %% moment variance in nAM
priormomvardefault=[1, 1, 1]; %% 



while adding_dips
    if dip_q>0, 
        msg_dip =['Add dipoles to ',num2str(dip_c),' or stop?'];
        dip_ch = 'Single|Symmetric Pair|Stop';
        dip_val = [1,2,0];
        def_opt=3;
    else
        msg_dip =['Add dipoles to model'];
        def_opt=1;
        dip_ch = 'Single|Symmetric Pair';
        dip_val = [1,2];
    end
    a_dip = spm_input(msg_dip,2+tr_q+dip_q,'b',dip_ch,dip_val,def_opt);
    if a_dip == 0
        adding_dips = 0;
    elseif a_dip == 1
    % add a single dipole to the model
        dip_q = dip_q+1;
        dip_pr(dip_q) = struct( 'a_dip',a_dip, ...
            'mu_w0',[],'mu_s0',[],'S_s0',eye(3),'S_w0',eye(3));
          %%  'ab20',[],'ab30',[]); %% ,'Tw', [],'Ts', []);
        % Location prior
        spr_q = spm_input('Location prior ?',1+tr_q+dip_q+1,'b', ...
                    'Informative|Non-info',[1,0],2);
        if spr_q
            % informative location prior
            str = 'Location prior';
            while 1
                s0mni = spm_input(str, 1+tr_q+dip_q+2,'e',[0 0 0])';
                
                %outside = ~ft_inside_headmodel(s0mni',mnivol);
                
                distcentre=sqrt(dot(s0mni'-mean(D.inv{D.val}.mesh.tess_mni.vert),s0mni'-mean(D.inv{D.val}.mesh.tess_mni.vert))); %% 
                outside=(distcentre>DISTTHRESH);
                    
                s0=D.inv{val}.datareg(modind).fromMNI*[s0mni' 1]';
                s0=s0(1:3);
                
                str2='Prior location variance (mm2)';
                diags_s0_mni = spm_input(str2, 1+tr_q+dip_q+2,'e',priorlocvardefault)';
                  
          
                S_s0_ctf=orM1*diag(diags_s0_mni)*orM1'; %% transform covariance
            
                
                %% need to leave diags(S0) free
              
                if all(~outside), break, end
                    str = 'Prior location must be inside head';
                  end
            dip_pr(dip_q).mu_s0 = s0;    
           else
            % no location  prior
            dip_pr(dip_q).mu_s0 = zeros(3,1);
            diags_s0_mni= nopriorlocvardefault';
            S_s0_ctf=diag(diags_s0_mni);
            end
        
        dip_pr(dip_q).S_s0=S_s0_ctf; %
        
        % Moment prior
        wpr_q = spm_input('Moment prior ?',1+tr_q+dip_q+spr_q+2,'b', ...
                    'Informative|Non-info',[1,0],2);
        if wpr_q
            % informative moment prior
            w0_mni= spm_input('Moment prior', ...
                                        1+tr_q+dip_q+spr_q+3,'e',[0 0 0])';
            str2='Prior moment variance (nAm2)';
            diags_w0_mni = spm_input(str2, 1+tr_q+dip_q+2,'e',priormomvardefault)';
            dip_pr(dip_q).mu_w0 =orM1*w0_mni;
            S_w0_ctf=orM1*diag(diags_w0_mni)*orM1';
            
             
        else
            % no location  prior
            dip_pr(dip_q).mu_w0 = zeros(3,1);
            S_w0_ctf= diag(nopriormomvardefault);
        end
        %% set up covariance matrix for orientation with no crosstalk terms (for single
        %% dip)
        dip_pr(dip_q).S_w0=S_w0_ctf;
        dip_c = dip_c+1;
    else
    % add a pair of symmetric dipoles to the model
        dip_q = dip_q+1;
        dip_pr(dip_q) = struct( 'a_dip',a_dip, ...
            'mu_w0',[],'mu_s0',[],'S_s0',eye(6),'S_w0',eye(6));
        %%...
          %  'ab20',[],'ab30',[]); %%,'Tw',eye(6),'Ts',eye(6));
        % Location prior
        spr_q = spm_input('Location prior ?',1+tr_q+dip_q+1,'b', ...
                    'Informative|Non-info',[1,0],2);
        if spr_q
            % informative location prior
             str = 'Location prior (one side only)';
            while 1
                s0mni = spm_input(str, 1+tr_q+dip_q+2,'e',[0 0 0])';
                syms0mni=s0mni;
                syms0mni(1)=-syms0mni(1);
                
                
                distcentre=sqrt(dot(s0mni'-mean(D.inv{D.val}.mesh.tess_mni.vert),s0mni'-mean(D.inv{D.val}.mesh.tess_mni.vert))); %% 
                outside=(distcentre>DISTTHRESH);
                
                s0=D.inv{val}.datareg(modind).fromMNI*[s0mni' 1]';
                
                s0sym=D.inv{val}.datareg(modind).fromMNI*[syms0mni' 1]';
                
                
                str2='Prior location variance (mm2)';
                tmp_diags_s0_mni = spm_input(str2, 1+tr_q+dip_q+2,'e',priorlocvardefault)';
                tmp_diags_s0_mni= [tmp_diags_s0_mni ; tmp_diags_s0_mni ];
                
              
                if all(~outside), break, end
                    str = 'Prior location must be inside head';
                  end
            
            dip_pr(dip_q).mu_s0 = [s0(1:3);s0sym(1:3)];    
        else                   
            % no location  prior
            dip_pr(dip_q).mu_s0 = zeros(6,1);
            tmp_diags_s0 = [nopriorlocvardefault';nopriorlocvardefault'];
            tmp_diags_s0_mni=[nopriorlocvardefault';nopriorlocvardefault'];
        end %% end of if informative prior
        %% setting up a covariance matrix where there is covariance between
        %% the x parameters negatively coupled, y,z positively.
         mni_dip_pr(dip_q).S_s0 = eye(length(tmp_diags_s0_mni)).*repmat(tmp_diags_s0_mni,1,length(tmp_diags_s0_mni));
         mni_dip_pr(dip_q).S_s0(4,1)=-mni_dip_pr(dip_q).S_s0(4,4); % reflect in x
         mni_dip_pr(dip_q).S_s0(5,2)=mni_dip_pr(dip_q).S_s0(5,5); % maintain y and z 
         mni_dip_pr(dip_q).S_s0(6,3)=mni_dip_pr(dip_q).S_s0(6,6);
         
         mni_dip_pr(dip_q).S_s0(1,4)=mni_dip_pr(dip_q).S_s0(4,1);
         mni_dip_pr(dip_q).S_s0(2,5)=mni_dip_pr(dip_q).S_s0(5,2);
         mni_dip_pr(dip_q).S_s0(3,6)=mni_dip_pr(dip_q).S_s0(6,3);
         %% transform to MEG space   
         
         %dip_pr(dip_q).S_s0(:,1:3)=orM1*mni_dip_pr(dip_q).S_s0(:,1:3)*orM1'; %% NEED TO LOOK AT THIS
         %dip_pr(dip_q).S_s0(:,4:6)=orM1*mni_dip_pr(dip_q).S_s0(:,4:6)*orM1';
         tmp1=orM1*mni_dip_pr(dip_q).S_s0(1:3,1:3)*orM1'; %% NEED TO LOOK AT THIS
         tmp2=orM1*mni_dip_pr(dip_q).S_s0(1:3,4:6)*orM1';
         tmp3=orM1*mni_dip_pr(dip_q).S_s0(4:6,4:6)*orM1';
         tmp4=orM1*mni_dip_pr(dip_q).S_s0(4:6,1:3)*orM1';
         dip_pr(dip_q).S_s0(1:3,1:3)=tmp1;
         dip_pr(dip_q).S_s0(1:3,4:6)=tmp2;
         dip_pr(dip_q).S_s0(4:6,4:6)=tmp3;
         dip_pr(dip_q).S_s0(4:6,1:3)=tmp4;
         
        % Moment prior
        wpr_q = spm_input('Moment prior ?',1+tr_q+dip_q+spr_q+2,'b', ...
                                           'Informative|Non-info',[1,0],2);
        if wpr_q
            % informative moment prior
            tmp= spm_input('Moment prior (right only)', ...
                                      1+tr_q+dip_q+spr_q+3,'e',[1 1 1])';
            tmp = [tmp ; tmp] ; tmp(4) = tmp(4);
            dip_pr(dip_q).mu_w0 = tmp;
            
            str2='Prior moment variance (nAm2)';
            diags_w0 = spm_input(str2, 1+tr_q+dip_q+spr_q+3,'e',priormomvardefault)';
            tmp_diags_w0=[diags_w0; diags_w0];
            
        else
            % no moment  prior
            dip_pr(dip_q).mu_w0 = zeros(6,1);
            tmp_diags_w0 = [nopriormomvardefault'; nopriormomvardefault'];
        end
        %dip_pr(dip_q).S_w0=eye(length(diags_w0)).*repmat(diags_w0,1,length(diags_w0));
        %% couple all orientations, except x, positively or leave for now...
          
                mni_dip_pr(dip_q).S_w0 = eye(length(tmp_diags_w0)).*repmat(tmp_diags_w0,1,length(tmp_diags_w0));
                
            mni_dip_pr(dip_q).S_w0(4,1)=-mni_dip_pr(dip_q).S_w0(4,4); % reflect x orientation
            mni_dip_pr(dip_q).S_w0(5,2)=mni_dip_pr(dip_q).S_w0(5,5); %
            mni_dip_pr(dip_q).S_w0(6,3)=mni_dip_pr(dip_q).S_w0(6,6); %
            mni_dip_pr(dip_q).S_w0(1,4)=-mni_dip_pr(dip_q).S_w0(4,1); %
            mni_dip_pr(dip_q).S_w0(2,5)=mni_dip_pr(dip_q).S_w0(5,2); %
            mni_dip_pr(dip_q).S_w0(3,6)=mni_dip_pr(dip_q).S_w0(6,3); %
            tmp1=orM1*mni_dip_pr(dip_q).S_w0(1:3,1:3)*orM1'; 
            tmp2=orM1*mni_dip_pr(dip_q).S_w0(1:3,4:6)*orM1';
            tmp3=orM1*mni_dip_pr(dip_q).S_w0(4:6,4:6)*orM1';
            tmp4=orM1*mni_dip_pr(dip_q).S_w0(4:6,1:3)*orM1';
            dip_pr(dip_q).S_w0(1:3,1:3)=tmp1;
            dip_pr(dip_q).S_w0(1:3,4:6)=tmp2;
            dip_pr(dip_q).S_w0(4:6,4:6)=tmp3;
            dip_pr(dip_q).S_w0(4:6,1:3)=tmp4;
         
        
        
        dip_c = dip_c+2;
    end
end

%str2='Data SNR (amp)';
% SNRamp = spm_input(str2, 1+tr_q+dip_q+2+1,'e',5)';
%SNRamp=3;
SNRamp=3;
 hE=log(SNRamp^2); %% expected log precision of data
 hC=1; % variability of the above precision
 
 
str2='Number of iterations';
 Niter = spm_input(str2, 1+tr_q+dip_q+2+2,'e',10)';
              
%%
% Get all the priors together and build structure to pass to inv_becd 

%============================

priors = struct('mu_w0',cat(1,dip_pr(:).mu_w0), ...
                'mu_s0',cat(1,dip_pr(:).mu_s0), ...
                'S_w0',blkdiag(dip_pr(:).S_w0),'S_s0',blkdiag(dip_pr(:).S_s0),'hE',hE,'hC',hC);
                
            
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
    
 %% set up figures 
    P.handles.hfig  = spm_figure('GetWin','Graphics');
    spm_clf(P.handles.hfig)
    P.handles.SPMdefaults.col = get(P.handles.hfig,'colormap');
    P.handles.SPMdefaults.renderer = get(P.handles.hfig,'renderer');
    set(P.handles.hfig,'userdata',P)
    dip_amp=[];
    for j=1:Niter,
        
     Pout(j) = spm_eeg_inv_vbecd(P);
     close(gcf);
     varresids(j)=var(Pout(j).y-Pout(j).ypost);
     pov(j)=100*(1-varresids(j)/var(Pout(j).y)); %% percent variance explained
     allF(j)=Pout(j).F;
     dip_mom=reshape(Pout(j).post_mu_w,3,length(Pout(j).post_mu_w)/3);
     dip_amp(j,:)=sqrt(dot(dip_mom,dip_mom));
    % display
     megloc=reshape(Pout(j).post_mu_s,3,length(Pout(j).post_mu_s)/3); % loc of dip (3 x n_dip)
     mniloc=D.inv{val}.datareg(modind).toMNI*[megloc;ones(1,size(megloc,2))]; %% actual MNI location (with scaling)
     megmom=reshape(Pout(j).post_mu_w,3,length(Pout(j).post_mu_w)/3); % moments of dip (3 x n_dip)
     megposvar=reshape(diag(Pout(j).post_S_s),3,length(Pout(j).post_mu_s)/3); %% estimate of positional uncertainty in three principal axes
     mnimom=orM1*megmom; %% convert moments into mni coordinates through a rotation (no scaling or translation)
     mniposvar=(orM1*sqrt(megposvar)).^2; %% convert pos variance into approx mni space by switching axes
     
     displayVBupdate2(Pout(j).y,pov,allF,Niter,dip_amp,mnimom,mniloc(1:3,:),mniposvar,P,j,[],Pout(j).F,Pout(j).ypost,[]);
     
    end; % for  j
    allF=[Pout.F];
    [maxFvals,maxind]=max(allF);


    P=Pout(maxind); %% take best F
    % Get the results out.
    inverse.pst = tb*1e3;
    inverse.F(ii) = P.F; % free energy
    
     megloc=reshape(P.post_mu_s,3,length(P.post_mu_s)/3); % loc of dip (3 x n_dip)
     meg_w=reshape(P.post_mu_w,3,length(P.post_mu_w)/3); % moments of dip (3 x n_dip)
     mni_w=orM1*meg_w; %% orientation in mni space
     mniloc=D.inv{val}.datareg(modind).toMNI*[megloc;ones(1,size(megloc,2))]; %% actual MNI location (with scaling)
     inverse.mniloc{ii}=mniloc(1:3,:);
    inverse.loc{ii} = megloc;
    
    
    
    
    inverse.j{ii} = P.post_mu_w; % dipole(s) orient/ampl, in 1 column in meg space
    inverse.jmni{ii} = reshape(mni_w,1,prod(size(mni_w)))'; % dipole(s) orient/ampl in mni space
    inverse.cov_loc{ii} = P.post_S_s; % cov matrix of source location
    inverse.cov_j{ii} = P.post_S_w; % cov matrix of source orient/ampl
    inverse.exitflag(ii) = 1; % Converged (1) or not (0)
    inverse.P{ii} = P; % save all kaboodle too.
    %% show final result
    pause(1);
    
    
    spm_clf(P.handles.hfig)
     
     megmom=reshape(Pout(maxind).post_mu_w,3,length(Pout(maxind).post_mu_w)/3); % moments of dip (3 x n_dip)
     %megposvar=reshape(diag(Pout(maxind).post_S_s),3,length(Pout(maxind).post_mu_s)/3); %% estimate of positional uncertainty in three principal axes
     mnimom=orM1*megmom; %% convert moments into mni coordinates through a rotation (no scaling or translation)
     
     
     longorM1=zeros(size(Pout(maxind).post_S_s,1));
     for k1=1:length(Pout(maxind).post_S_s)/3;
        longorM1((k1-1)*3+1:k1*3,(k1-1)*3+1:k1*3)=orM1;
     end; % for k1
     S0_mni=longorM1*Pout(maxind).post_S_s*longorM1';
     mniposvar=diag(S0_mni); %% convert pos variance into approx mni space by switching axes
     mniposvar=reshape(mniposvar,3,length(Pout(maxind).post_S_s)/3);
     
     displayVBupdate2(Pout(maxind).y,pov,allF,Niter,dip_amp,mnimom,mniloc(1:3,:),mniposvar,P,j,[],Pout(maxind).F,Pout(maxind).ypost,maxind);
     
    %displayVBupdate2(Pout(maxind).y,pov,allF,Niter,dip_amp,mniloc,Pout(maxind).post_mu_s,Pout(maxind).post_S_s,P,j,[],Pout(maxind).F,Pout(maxind).ypost,maxind,D);
  % 
end

%  [tmp] = ft_compute_leadfield([-42 15 48], P.forward.sens, P.forward.vol,'reducerank',2);
%  estfield=tmp*[-21 -8 -45]';
%  plot(P.y,estfield./1000,'x')
%  axis equal
 
D.inv{val}.inverse = inverse;

%%
% Save results and display
%-------------------------
save(D)


return





function [P] = displayVBupdate2(y,pov_iter,F_iter,maxit,dipamp_iter,mu_w,mu_s,diagS_s,P,it,flag,F,yHat,maxind)


%% yHat is estimate of y based on dipole position
if ~exist('flag','var')
    flag = [];
end
if ~exist('maxind','var')
    maxind = [];
end

if isempty(flag) || isequal(flag,'ecd')
    % plot dipoles
    try
        opt.ParentAxes = P.handles.axesECD;
        opt.hfig = P.handles.hfig;
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
            'visible','off',...
            'deleteFcn',@back2defaults);
        opt.ParentAxes = P.handles.axesECD;
        opt.hfig = P.handles.hfig;
    end
    w = reshape(mu_w,3,[]);
    s = reshape(mu_s, 3, []);
%     mesh.faces=D.inv{D.val}.forward.mesh.face; %% in ctf space
%     mesh.vertices=D.inv{D.val}.forward.mesh.vert; %% in ctf space
%     [out] = spm_eeg_displayECD_ctf(...
%         s,w,reshape(diag(S_s),3,[]),[],mesh,opt);
       [out] = spm_eeg_displayECD(...
         s,w,diagS_s,[],opt);
    
        P.handles.hp = out.handles.hp;
        P.handles.hq = out.handles.hq;
        P.handles.hs = out.handles.hs;
        P.handles.ht = out.handles.ht;
     
end

% plot data and predicted data
pos = P.forward.sens.prj;
ChanLabel = P.channels;
in.f = P.handles.hfig;
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
if isempty(flag) || isequal(flag,'data') || isequal(flag,'ecd')
    %yHat = P.gmn*mu_w;
    miY = min([yHat;y]);
    maY = max([yHat;y]);
    try
        P.handles.axesYhat;
        d = get(P.handles.axesYhat,'userdata');
        yHat = yHat(d.goodChannels);
        clim = [min(yHat(:))-( max(yHat(:))-min(yHat(:)) )/63,...
            max(yHat(:))];
        ZI = griddata(...
            d.interp.pos(1,:),d.interp.pos(2,:),full(double(yHat)),...
            d.interp.XI,d.interp.YI);
        set(d.hi,'Cdata',flipud(ZI));
        caxis(P.handles.axesYhat,clim);
        delete(d.hc)
        [C,d.hc] = contour(P.handles.axesYhat,flipud(ZI),...
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

end


if isempty(flag) || isequal(flag,'var')
    % plot precision hyperparameters
    try
        P.handles.axesVar1;
    catch
        figure(P.handles.hfig)
        P.handles.axesVar1 = axes(...
            'Position',[0.05 0.05 0.25 0.2],...
            'NextPlot','replace',...
            'box','on');
    end
    plot(P.handles.axesVar1,F_iter,'o-');
    if ~isempty(maxind),
        hold on;
        h=plot(P.handles.axesVar1,maxind,F_iter(maxind),'rd');
        set(h,'linewidth',4);
        end;
    set(P.handles.axesVar1,'Xlimmode','manual');
    set(P.handles.axesVar1,'Xlim',[1 maxit]);
    set(P.handles.axesVar1,'Xtick',1:maxit);
    set(P.handles.axesVar1,'Xticklabel',num2str([1:maxit]'));
    set(P.handles.axesVar1,'Yticklabel','');
    title(P.handles.axesVar1,'Free energy ')
    axis(P.handles.axesVar1,'square');
    set(P.handles.axesVar1,'Ylimmode','auto'); %,'tight')
    
    grid(P.handles.axesVar1,'on')

    try
        P.handles.axesVar2;
    catch
        figure(P.handles.hfig)
        P.handles.axesVar2 = axes(...
            'Position',[0.37 0.05 0.25 0.2],...
            'NextPlot','replace',...
            'box','on');
    end
    
    
    plot(P.handles.axesVar2,pov_iter,'*-') 
    if ~isempty(maxind),
        hold on;
        h=plot(P.handles.axesVar2,maxind,pov_iter(maxind),'rd');
        set(h,'linewidth',4);
        end;   
    set(P.handles.axesVar2,'Xlimmode','manual');
    set(P.handles.axesVar2,'Xlim',[1 maxit]);
    set(P.handles.axesVar2,'Xtick',1:maxit);
    set(P.handles.axesVar2,'Xticklabel',num2str([1:maxit]'));
    set(P.handles.axesVar2,'Ylimmode','manual'); %,'tight')
     set(P.handles.axesVar2,'Ylim',[0 100]);
    set(P.handles.axesVar2,'Ytick',[0:20:100]);
    set(P.handles.axesVar2,'Yticklabel',num2str([0:20:100]'));
 
    %set(P.handles.axesVar2,'Yticklabel','');
    title(P.handles.axesVar2,'Percent variance explained');
    axis(P.handles.axesVar2,'square');
    
    
    grid(P.handles.axesVar2,'on')
    

    try
        P.handles.axesVar3;
    catch
        figure(P.handles.hfig)
        P.handles.axesVar3 = axes(...
            'Position',[0.72 0.05 0.25 0.2],...
            'NextPlot','replace',...
            'box','on');
    end
    plot(P.handles.axesVar3,1:it,dipamp_iter','o-');
      if ~isempty(maxind),
        hold on;
        h=plot(P.handles.axesVar3,maxind,dipamp_iter(maxind,:)','rd');
        set(h,'linewidth',4);
        end;
    
    set(P.handles.axesVar3,'Xlimmode','manual');
    set(P.handles.axesVar3,'Xlim',[1 maxit]);
    set(P.handles.axesVar3,'Xtick',1:maxit);
    set(P.handles.axesVar3,'Xticklabel',num2str([1:maxit]'));
    set(P.handles.axesVar3,'Yticklabel','');
    title(P.handles.axesVar3,'Dipole amp (nAm) ')
    axis(P.handles.axesVar3,'square');
    set(P.handles.axesVar3,'Ylimmode','auto'); %,'tight')
    grid(P.handles.axesVar3,'on')

    
    
end



if ~isempty(flag) && (isequal(flag,'ecd') || isequal(flag,'mGN') )
    try
        P.handles.hte(2);
    catch
        figure(P.handles.hfig)
        P.handles.hte(2) = uicontrol('style','text',...
            'units','normalized',...
            'position',[0.2,0.91,0.6,0.02],...
            'backgroundcolor',[1,1,1]);
    end
    set(P.handles.hte(2),'string',...
        ['ECD locations: Modified Gauss-Newton scheme... ',num2str(floor(P.pc)),'%'])
else
    try
        set(P.handles.hte(2),'string','VB updates on hyperparameters')
    end       
end

try
    P.handles.hte(1);
catch
    figure(P.handles.hfig)
    P.handles.hte(1) = uicontrol('style','text',...
        'units','normalized',...
        'position',[0.2,0.94,0.6,0.02],...
        'backgroundcolor',[1,1,1]);
end
try
    set(P.handles.hte(1),'string',...
        ['Model evidence: p(y|m) >= ',num2str(F(end),'%10.3e\n')])
end

try
    P.handles.hti;
catch
    figure(P.handles.hfig)
    P.handles.hti = uicontrol('style','text',...
        'units','normalized',...
        'position',[0.3,0.97,0.4,0.02],...
        'backgroundcolor',[1,1,1],...
        'string',['VB ECD inversion: trial #',num2str(P.ltr(P.ii))]);
end
drawnow

function back2defaults(e1,e2)
hf = spm_figure('FindWin','Graphics');
P = get(hf,'userdata');
try
    set(hf,'colormap',P.handles.SPMdefaults.col);
    set(hf,'renderer',P.handles.SPMdefaults.renderer);
end



