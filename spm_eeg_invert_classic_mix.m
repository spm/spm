function [D] = spm_eeg_invert_classic_mix(D,val,Qpriors,surfind,ugainfiles)
%function [D] = spm_eeg_invert_classic_mix(D,val,Qpriors,surfind,ugainfiles)
%%
% ReML inversion of multiple posterior current variances from previous
% iterations spm_eeg_invert_classic or spm_eeg_invert
% ReML estimation of regularisation hyperparameters using the
% spatiotemporal hierarchy implicit in EEG/MEG data
%
% D contains the data and inversion parameters (see
% spm_eeg_invert.m/spm_eeg_invert_classic.m)
% val the inversion index 
% Qpriors is N solutions of rows by Nd variance estimates
% surfind is N solutions long and contains indices into ugainfiles to these priors with different lead field structures
%% ugainfiles are the SPMgain matrices for the different surfaces
%
% Output D will have a solution which is optimal REML mixture of Qpriors
%
% Created by:   
%               Gareth Barnes - g.barnes@ucl.ac.uk
%
% This version is for single subject single modality analysis and therefore
% contains none of the associated scaling factors.
%
% $Id: spm_eeg_invert_classic_mix.m 6077 2014-06-30 16:55:03Z spm $



Nl = length(D);



if Nl>1,
    error('function only defined for a single subject');
end;

% D - SPM data structure
%==========================================================================
if nargin > 1
    D.val = val;
elseif ~isfield(D, 'val')
    D.val = 1;
end


val=D.val;

inverse   = D.inv{val}.inverse;

% forward model
%--------------------------------------------------------------------------


% defaults
%--------------------------------------------------------------------------

try, type = inverse.type;   catch, type = 'GS';     end
try, s    = inverse.smooth; catch, s    = 0.6;      end
try, Np   = inverse.Np;     catch, Np   = 256;      end
try, Nr   = inverse.Nr;     catch, Nr   = 16;       end %% requested number of temporal modes, could be changed depending on svd
try, xyz  = inverse.xyz;    catch, xyz  = [0 0 0];  end
try, rad  = inverse.rad;    catch, rad  = 128;      end
try, hpf  = inverse.hpf;    catch, hpf  = 48;       end %% need to one day put these the correct way round
try, lpf  = inverse.lpf;    catch, lpf  = 0;        end
try, sdv  = inverse.sdv;    catch, sdv  = 4;        end
try, Han  = inverse.Han;    catch, Han  = 1;        end
try, woi  = inverse.woi;    catch, woi  = [];       end
try, Nm   = inverse.Nm;     catch, Nm   = [];       end
try, Nt   = inverse.Nt;     catch, Nt   = [];       end %% fixed number of temporal modes
try, Ip   = inverse.Ip;     catch, Ip   = [];       end
try, SHUFFLELEADS=inverse.SHUFFLELEADS;catch, SHUFFLELEADS=[];end




% defaults
%--------------------------------------------------------------------------
%type = inverse.type;    % Type of inversion scheme


% get specified modalities to invert (default to all)

%--------------------------------------------------------------------------
modalities = D.inv{val}.forward.modality;       % MEG in this case

Nmax  = 16;         % max number of temporal modes

% check lead fields and get number of dipoles (Nd) and channels (Nc)
%==========================================================================

if inverse.Nm~=length(inverse.Ic{1}),
    disp('Using reduced spatial modes');
    if length(unique(surfind))>1,
        error('Cannot merge different surface files with reduced spatial modes at the moment');
    end;
    U=inverse.U{1};
else
    disp('Using all spatial modes');
    U=eye(inverse.Nm);
end;

    

fprintf('Checking leadfields')



usurfind=unique(surfind);
Nsurf=length(usurfind);
if Nsurf~=size(ugainfiles,1); %% number of different surfaces
    error('need a lead field file for each surface');
end;

L=[];
g_ind=zeros(Nsurf,2); %% index start and end of each set of lead fields
for f=1:Nsurf; %% load in  all leadfield
    
    fname = deblank(ugainfiles(f,:));
    G = load(fullfile(fname)); % Relative path
    
    
    label = G.label;
    g_ind(f,1)=size(L,2)+1; %% number of lead fields per surface
    L     =[L G.G];
    g_ind(f,2)=size(L,2); %% number of lead fields per surface
    clear G;
    
end;




if size(modalities,1)>1,
    error('not defined for multiple modalities');
end;
Ic  = setdiff(D.indchantype(modalities), badchannels(D));
Nd    = size(L,2);      % Number of dipoles

fprintf(' - done\n')



% check for (e.g., empty-room) sensor components (in Qe)
%==========================================================================
QE = 1;                     % No empty room noise measurement


%==========================================================================
% Spatial projectors - keep these to map directly to channels for now
%==========================================================================



A=U;
UL=U*L;

clear L

Nm    = size(UL,1);         % Number of spatial projectors


% Report
%----------------------------------------------------------------------
fprintf('Using %d spatial modes',Nm)

% No dipole is eliminated
%--------------------------------------------------------------------------
Is    = 1:Nd;               % Accepted dipoles
%Ns    = length(Is);         % Ns = Nd in this case


%==========================================================================
% Temporal projector
%==========================================================================
AY    = {};                                      % pooled response for MVB
AYYA  = 0;                                       % pooled response for ReML

% Time-window of interest
%----------------------------------------------------------------------

if isempty(woi)
    w      = 1000*[min(D.time) max(D.time)];
else
    w=woi; %% in milliseconds
end;

It     = (w/1000 - D.timeonset)*D.fsample + 1;
It     = max(1,It(1)):min(It(end), length(D.time));
It     = fix(It);

% Peristimulus time
%----------------------------------------------------------------------
pst    = 1000*D.time;                   % peristimulus time (ms)
pst    = pst(It);                       % windowed time (ms)
dur    = (pst(end) - pst(1))/1000;      % duration (s)
dct    = (It - It(1))/2/dur;            % DCT frequencies (Hz)
Nb     = length(It);                    % number of time bins

% Serial correlations
%----------------------------------------------------------------------
K      = exp(-(pst - pst(1)).^2/(2*sdv^2));
K      = toeplitz(K);
qV     = sparse(K*K');

% Confounds and temporal subspace
%----------------------------------------------------------------------

T      = spm_dctmtx(Nb,Nb);         % use plot(T) here!

j      = find( (dct >= lpf) & (dct <= hpf) ); %% THis is the wrong way round but leave for nowfor compatibility with spm_eeg_invert
T      = T(:,j);                    % Apply the filter to discrete cosines
dct    = dct(j);                    % Frequencies accepted

%% Hanning window
%----------------------------------------------------------------------

if Han
    W  = sparse(1:Nb,1:Nb,spm_hanning(Nb)); %% use hanning unless specified
else
    W=1;
end;




% get trials or conditions
%----------------------------------------------------------------------
try
    trial = D.inv{D.val}.inverse.trials;
catch
    trial = D.condlist;
end
Ntrialtypes=length(trial);
% get temporal covariance (Y'*Y) to find temporal modes
%======================================================================
%MY    = cell(Nmod,1);                        % mean response
YTY   = sparse(0);                           % accumulator


% get (spatially aligned) data
%------------------------------------------------------------------

YY    = 0;
%    MY{m} = 0;
N=0;
badtrialind=D.badtrials;
for j = 1:Ntrialtypes,                          % pool over conditions
    c     = D.indtrial(trial{j});     % and trials
    c=setxor(c,badtrialind);
    length(c)
    Nk    = length(c);
    for k = 1:Nk
        Y     = A*D(Ic,It,c(k));
        
        YY    = YY + Y'*Y;
        N     = N + 1;
    end
end
YY=YY./N;




% Apply any Hanning and filtering
%------------------------------------------------------------------
YY         = W'*YY*W;     % Hanning
YTY         = T'*YY*T;     % Filter


%======================================================================

if isempty(Nt),
    
    [U E]  = spm_svd(YTY,exp(-8));          % get temporal modes
    if isempty(U),
        warning('nothing found using spm svd, using svd');
        [U E]  = svd(YTY);          % get temporal modes
    end;
    E      = diag(E)/trace(YTY);            % normalise variance
    Nr     = min(length(E),Nmax);           % number of temporal modes
    Nr=max(Nr,1); %% use at least one mode
else
    [U E]  = svd(YTY);          % get temporal modes
    
    E      = diag(E)/trace(YTY);            % normalise variance
    disp('Fixed number of temporal modes');
    Nr=Nt;
end;

V      = U(:,1:Nr);                     % temporal modes
VE     = sum(E(1:Nr));                  % variance explained

fprintf('Using %i temporal modes, ',Nr)
fprintf('accounting for %0.2f percent average variance\n',full(100*VE))

% projection and whitening
%----------------------------------------------------------------------
S      = T*V;                           % temporal projector
Vq     = S*pinv(S'*qV*S)*S';            % temporal precision


% get spatial covariance (Y*Y') for Gaussian process model
%======================================================================

%======================================================================
%======================================================================

% loop over Ntrialtypes trial types
%----------------------------------------------------------------------
UYYU = 0;
AYYA=0;
Nn    =0;                             % number of samples
AY={};
Ntrials=0;
for j = 1:Ntrialtypes,
    
    UY{j} = sparse(0);
    c       = D.indtrial(trial{j});
    Nk      = length(c);
    
    % loop over epochs
    %------------------------------------------------------------------
    for k = 1:Nk
        
        % stack (scaled aligned data) over modalities
        %--------------------------------------------------------------
        
        Y       = D(Ic,It,c(k))*S;
        Y=A*Y; %% CHANGE
        
        
        
        % accumulate first & second-order responses
        %--------------------------------------------------------------
        Nn       = Nn + Nr;         % number of samples
        %Y           = spm_cat(MY);           % contribution to ERP, Y and MY are the same for 1 modality
        YY          = Y*Y';                  % and covariance
        Ntrials=Ntrials+1;
        %YYep{k}=YY;%% one for each trial and condition
        
        % accumulate statistics (subject-specific)
        %--------------------------------------------------------------
        UY{j}     = UY{j} + Y;           % condition-specific ERP
        UYYU     = UYYU + YY;          % subject-specific covariance
        
        % and pool for optimisation of spatial priors over subjects
        %--------------------------------------------------------------
        AY{end + 1} = Y;                     % pooled response for MVB
        AYYA        = AYYA    + YY;          % pooled response for ReML
        
    end
end

AYYA=AYYA./Nn;
AY=spm_cat(AY);


% assuming equal noise over subjects (Qe) and modalities AQ
%--------------------------------------------------------------------------
AQeA   = A*QE*A';           % Note that here it is A*A'
Qe{1}  = AQeA/(trace(AQeA)); % it means IID noise in virtual sensor space


% Inverse solution
%==========================================================================
QP     = {};
LQP    = {};
LQPL   = {};


for j=1:size(Qpriors,1),
    % Accumulate empirical priors (New set of patches for the second inversion)
    %------------------------------------------------------------------
    
    disp(sprintf('Adding prior set %d of %d',j,size(Qpriors,1)));
    dum=sparse(zeros(1,Nd));
    
    dum(g_ind(surfind(j),1):g_ind(surfind(j),2))=Qpriors(j,:);
    
    qp=sparse(diag(dum)); %% qp is Nd*Nd
    QP{j}=sparse(diag(qp));
    LQP{j}=UL*qp;
    LQPL{j} = LQP{j}*UL';
    clear qp
end;

%==========================================================================
% Step 2: Re-estimate for each subject separately (fusing all modalities)
%==========================================================================

fprintf('Inverting subject 1\n')

% generate sensor component (Qe) per modality
%----------------------------------------------------------------------
AQeA  = A*QE*A';                % Again it is A*A'
AQ    = AQeA/(trace(AQeA));


% using spatial priors
%----------------------------------------------------------------------
Np    = length(LQPL);       % Final number of priors
Ne    = length(Qe);         % Sensor noise prior
Q     = [Qe LQPL];


% re-do ReML (with informative hyperpriors)
% Here is performed the second inversion
%======================================================================

Q0          = exp(-2)*trace(AYYA)*AQ/trace(AQ);
[Cy,h,Ph,F] = spm_reml_sc(AYYA,[],Q,1,-4,16,Q0);


% Data ID
%----------------------------------------------------------------------
% When should I use the data ID?
% For comparison purposes it is necessary to guarantee that the IDs have
% the same value.
%
% When using the same dataset but different lead field matrices is a good
% example of fail, because the spatial projector will generate different
% virtual sensors, and therefore different data for the inversion.
ID    = spm_data_id(YY);

% Covariance: sensor space - Ce and source space - L*Cp
%----------------------------------------------------------------------
Cp    = sparse(0);
LCp   = sparse(0);
hp    = h(Ne + (1:Np));
for j = 1:Np
    Cp  =  Cp + hp(j)*QP{j};
    LCp = LCp + hp(j)*LQP{j};
end

% MAP estimates of instantaneous sources
%======================================================================
% This is equivalent to M = Cp*UL'*inv(Qe + UL*Cp*UL'))
% with Cp the posterior source covariance (with optimal h values)
M     = LCp'/Cy;

% conditional variance (leading diagonal)
% Cq    = Cp - Cp*L'*iC*L*Cp;
%----------------------------------------------------------------------
Cq    = Cp - sum(LCp.*M')';

% evaluate conditional expectation
%----------------------------------------------------------------------
% evaluate conditional expectation (of the sum over trials)
%----------------------------------------------------------------------
SSR   = 0;
SST   = 0;
J     = {};

for j = 1:Ntrialtypes
    
    % trial-type specific source reconstruction
    %------------------------------------------------------------------
    J{j} = M*UY{j};
    
    % sum of squares
    %------------------------------------------------------------------
    SSR  = SSR + sum(var((UY{j} - UL*J{j}))); %% changed variance calculation
    SST  = SST + sum(var( UY{j}));
    
end

% J = M*Y;
%
% % sum of squares
% %------------------------------------------------------------------
% SSR  = sum(var((Y - UL*J),0,2));
% SST  = sum(var(Y,0,2));

% accuracy; signal to noise (over sources)
%======================================================================
R2   = 100*(SST - SSR)/SST;
fprintf('Percent variance explained %.2f (%.2f)\n',full(R2),full(R2*VE));

% Save results
% DEMO: WARNING! These results are not coincident in format with
%                those generated in the SPM8
%======================================================================
inverse.type   = type;                 % inverse model
inverse.M      = M;                    % MAP projector (reduced)
inverse.J   = J;                    % Conditional expectation
inverse.Y      = Y;                    % ERP data (reduced)
inverse.L      = UL;                   % Lead-field (reduced)
inverse.qC     = Cq;                   % spatial covariance
inverse.qV     = Vq;                   % temporal correlations
inverse.T      = S;                    % temporal projector
inverse.U      = {A};                    % spatial projector
inverse.Is     = Is;                   % Indices of active dipoles
inverse.It     = It;                   % Indices of time bins
try
    inverse.Ic{1}     = Ic;                   % Indices of good channels
catch
    inverse.Ic    = Ic;                   % Indices of good channels
end;
inverse.Nd     = Nd;                   % number of dipoles
inverse.pst    = pst;                  % peristimulus time
inverse.dct    = dct;                  % frequency range
inverse.F      = F;                    % log-evidence
inverse.ID     = ID;                   % data ID
inverse.R2     = R2;                   % variance explained (reduced)
inverse.VE     = R2*VE;                % variance explained
inverse.woi    = w;                    % time-window inverted

inverse.modality = modalities;         % modalities inverted


% save in struct
%----------------------------------------------------------------------
D.inv{val}.inverse = inverse;
D.inv{val}.method  = 'Imaging';

% display
%======================================================================
spm_eeg_invert_display(D);
drawnow

return
