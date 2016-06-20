function [D] = spm_eeg_invert_classic(D,val)
%
%% A trimmed down version of the spm_eeg_invert() routine
%% this version only handles single subject single modality data
%% the removal of many scaling factors makes it easier to compare between forward models
% ReML inversion of multiple forward models for EEG-MEG
% FORMAT [D] = spm_eeg_invert_classic(D)
% ReML estimation of regularisation hyperparameters using the
% spatiotemporal hierarchy implicit in EEG/MEG data
%
% Requires:
% D{i}.inv{val}.inverse:
%
%     inverse.modality - modality to use in case of multimodal datasets
%
%     inverse.trials - D.events.types to invert
%     inverse.type   - 'GS' Greedy search on MSPs
%                      'ARD' ARD search on MSPs
%                      'MSP' GS and ARD multiple sparse priors
%                      'LOR' LORETA-like model
%                      'IID' minimum norm
%                       'EBB' for empirical bayes beamformer
%     inverse.woi    - time window of interest ([start stop] in ms)
%     inverse.lpf    - band-pass filter - low frequency cut-off (Hz)
%     inverse.hpf    - band-pass filter - high frequency cut-off (Hz)
%     inverse.Han    - switch for Hanning window
%     inverse.xyz    - (n x 3) locations of spherical VOIs
%     inverse.rad    - radius (mm) of VOIs
%
%     inverse.Nm     - maximum number of channel modes
%     inverse.Nr     - maximum number of temporal modes
%     inverse.Np     - number of sparse priors per hemisphere
%     inverse.smooth - smoothness of source priors (0 to 1)
%     inverse.Na     - number of most energetic dipoles
%     inverse.sdv    - standard deviations of Gaussian temporal correlation
%     inverse.Qe     - any sensor error components (e.g. empty-room data)
%     inverse.Qe0     - minimum amount of sensor noise power relative to
%                        signal eg 0.1 would correspond to power SNR of 10.0
%     inverse.A       - predefined spatial modes (Nchans*Nmodes) to project
%                       sensor data through
%
% Evaluates:
%
%     inverse.M      - MAP projector (reduced)
%     inverse.J{i}   - Conditional expectation (i conditions) J = M*U*Y
%     inverse.L      - Lead field (reduced UL := U*L)
%     inverse.qC     - spatial covariance
%     inverse.qV     - temporal correlations
%     inverse.T      - temporal projector
%     inverse.U(j)   - spatial projector (j modalities) - derived from data
%     inverse.A      - pre-specified spatial projector
%     inverse.Y{i}   - reduced data (i conditions) UY = UL*J + UE
%     inverse.Is     - Indices of active dipoles
%     inverse.It     - Indices of time bins
%     inverse.Ic{j}  - Indices of good channels (j modalities)
%     inverse.Nd     - number of dipoles
%     inverse.pst    - peristimulus time
%     inverse.dct    - frequency range
%     inverse.F      - log-evidence
%     inverse.VE     - variance explained in spatial/temporal subspaces (%)
%     inverse.R2     - variance in subspaces accounted for by model (%)
%     inverse.scale  - scaling of data for each of j modalities
%__________________________________________________________________________
%
% Created by:   Jose David Lopez - ralph82co@gmail.com
%               Gareth Barnes - g.barnes@fil.ion.ucl.ac.uk
%               Vladimir Litvak - litvak.vladimir@gmail.com
%
%
% This version is for single subject single modality analysis and therefore
% contains none of the associated scaling factors.
% No symmetric priors are used in this implementation (just single patches)
% There is an option for a Beamforming prior : inversion type 'EBB'
% also added new beamforming method- using GS rather than ARD- from Juan David Martinez Vargas 'EBBgs'

%%The code was used in
%% L�pez, J. D., Penny, W. D., Espinosa, J. J., Barnes, G. R. (2012).
% A general Bayesian treatment for MEG source reconstruction incorporating lead field uncertainty.
% Neuroimage 60(2), 1194-1204 doi:10.1016/j.neuroimage.2012.01.077.

% $Id: spm_eeg_invert_classic.m 6815 2016-06-20 09:08:18Z gareth $



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
try, QE    = inverse.QE;     catch,  QE=1;          end         %  empty room noise measurement
try, Qe0   = inverse.Qe0;     catch, Qe0   = exp(-5);       end  %% set noise floor at 1/100th signal power i.e. assume amplitude SNR of 10
try, inverse.A;     catch, inverse.A   = [];       end %% orthogonal channel modes

try, SHUFFLELEADS=inverse.SHUFFLELEADS;catch, SHUFFLELEADS=0;end; %% ONLY FOR TESTING - destroyes correspondence between geometry and data


% defaults
%--------------------------------------------------------------------------
type = inverse.type;    % Type of inversion scheme


% get specified modalities to invert (default to all)

%--------------------------------------------------------------------------
modalities = D.inv{val}.forward.modality;       % MEG in this case

Nmax  = 16;         % max number of temporal modes

% check lead fields and get number of dipoles (Nd) and channels (Nc)
%==========================================================================

fprintf('Checking leadfields')
[L,D] = spm_eeg_lgainmat(D);    % Generate/load lead field
Nd=size(L,2);

if ~isempty(Ip)
    Np   = length(Ip);              % Number of priors/3 for GS, ARD, MSP
else
    Ip=ceil([1:Np]*Nd/Np);
end;

persistent permind;

if SHUFFLELEADS,
    if isempty(permind),
        permind=randperm(size(L,1));
    end;
    L=L(permind,:);
    warning('PERMUTING LEAD FIELDS !');
    permind(1:3)
end;

% Check gain or lead-field matrices
%------------------------------------------------------------------

if size(modalities,1)>1,
    error('not defined for multiple modalities');
end;
Ic  = setdiff(D.indchantype(modalities), badchannels(D));
Nd    = size(L,2);      % Number of dipoles

fprintf(' - done\n')

if s>=1,
    smoothtype='mesh_smooth',
else
    smoothtype='msp_smooth'
end;
vert  = D.inv{val}.mesh.tess_mni.vert;
face  = D.inv{val}.mesh.tess_mni.face;
M1.faces=face;
M1.vertices=vert;

switch smoothtype,
    case 'mesh_smooth',
        fprintf('Using SPM smoothing for priors:')
        %--------------------------------------------------------------------------
        
        Qi    = speye(Nd,Nd);
        [QG]=spm_mesh_smooth(M1,Qi,round(s));
        QG    = QG.*(QG > exp(-8));     % Eliminate small values
        
        QG    = QG*QG;              % Guarantee positive semidefinite matrix
        disp('Normalising smoother');
        QG=QG./repmat(sum(QG,2),1,size(QG,1));
        
        
        
    case 'msp_smooth'
        fprintf('Computing Green function from graph Laplacian to smooth priors:')
        %--------------------------------------------------------------------------
        
        A     = spm_mesh_distmtx(struct('vertices',vert,'faces',face),0);
        
        GL    = A - spdiags(sum(A,2),0,Nd,Nd);
        GL    = GL*s/2;
        Qi    = speye(Nd,Nd);
        QG    = sparse(Nd,Nd);
        
        for i = 1:8
            QG = QG + Qi;
            Qi = Qi*GL/i;
        end
        
        QG    = QG.*(QG > exp(-8));
        QG    = QG*QG;
end;


clear Qi A GL
fprintf(' - done\n')




% check for (e.g., empty-room) sensor components (in Qe)
%==========================================================================


%==========================================================================
% Spatial projectors (adjusting for different Lead-fields)
%==========================================================================

fprintf('Optimising and aligning spatial modes ...\n')

% eliminate low SNR spatial modes
%------------------------------------------------------------------

if isempty(inverse.A), % no spatial modes prespecified
    if isempty(Nm), %% number of modes not specifiedd
        [U,ss,vv]    = spm_svd((L*L'),exp(-16));
        A     = U';                 % spatial projector A
        UL    = A*L;
        
    else % number of modes pre-specified
        [U,ss,vv]    = spm_svd((L*L'),0);
        if length(ss)<Nm,
            disp('number available');
            length(ss)
            error('Not this many spatial modes in lead fields');
            
        end;
        
        ss=ss(1:Nm);
        disp('using preselected number spatial modes !');
        A     = U(:,1:Nm)';                 % spatial projector A
        UL    = A*L;
    end;
else %% U was specified in input
    disp('Using pre-specified spatial modes');
    if isempty(Nm),
        error('Need to specify number of spatial modes if U is prespecified');
    end;
    %
    A=inverse.A;
    UL=A*L;
end;

Nm    = size(UL,1);         % Number of spatial projectors

clear ss vv

% Report
%----------------------------------------------------------------------
fprintf('Using %d spatial modes',Nm)

% None dipole is eliminated
%--------------------------------------------------------------------------
Is    = 1:Nd;               % Accepted dipoles
Ns    = length(Is);         % Ns = Nd in this case


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
disp(sprintf('Number of samples %d',length(It)));

    

% Peristimulus time
%----------------------------------------------------------------------
pst    = 1000*D.time;                   % peristimulus time (ms)
pst    = pst(It);                       % windowed time (ms)
dur    = (pst(end) - pst(1))/1000;      % duration (s)
dct    = (It - It(1))/2/dur;            % DCT frequencies (Hz)
Nb     = length(It);                    % number of time bins

% Serial correlations
%----------------------------------------------------------------------
K      = exp(-(pst - pst(1)).^2/(2*sdv^2)); %% sdv set to 4 by default
K      = toeplitz(K);
qV     = sparse(K*K'); %% Samples* samples covariance matrix- assumes smooth iid

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

YY    = 0;

N=0;

badtrialind=D.badtrials;
Ik=[]; %% keep a record of trials used
for j = 1:Ntrialtypes,                          % pool over conditions
    c     = D.indtrial(trial{j});     % and trials
    [c1,ib]=intersect(c,badtrialind); %% remove bad trials ib if there are any
    c=c(setxor(1:length(c),ib));
    Ik=[Ik c];
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

if isempty(Nt), %% automatically assign appropriate number of temporal modes    
    [U E]  = spm_svd(YTY,exp(-8));          % get temporal modes
    if isempty(U), %% fallback
        warning('nothing found using spm svd, using svd');
        [U E]  = svd(YTY);          % get temporal modes
    end;
    E      = diag(E)/trace(YTY);            % normalise variance
    Nr     = min(length(E),Nmax);           % number of temporal modes
    Nr=max(Nr,1); %% use at least one mode
else %% use predefined number of modes
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
    [c1,ib]=intersect(c,badtrialind); %% remove bad trials ib if there are any
    c=c(setxor(1:length(c),ib));
    Nk    = length(c);
       
    % loop over epochs
    %------------------------------------------------------------------
    for k = 1:Nk
        
        % stack (scaled aligned data) over modalities
        %--------------------------------------------------------------
        
        Y       = D(Ic,It,c(k))*S; %% in temporal subspace
        Y=A*Y; %%  in spatial subspace
  
        
        % accumulate first & second-order responses
        %--------------------------------------------------------------
        Nn       = Nn + Nr;         % number of samples
        
        YY          = Y*Y';                  % and covariance
        Ntrials=Ntrials+1;
        
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

AY=spm_cat(AY); %% goes to MVB/GS algorithm

ID    = spm_data_id(AY); %% get a unique ID for these filtered data

% assuming equal noise over subjects (Qe) and modalities AQ
%--------------------------------------------------------------------------
AQeA   = A*QE*A';           % Note that here it is A*A'
Qe{1}  = AQeA/(trace(AQeA)); % it means IID noise in virtual sensor space

%Q0          = Qe0*trace(AYYA)*Qe{1}*Nr; %% fixed (min) level of sensor space variance- this is divided by Nr later in spm_reml_sc
Q0          = Qe0*trace(AYYA)*Qe{1}./sum(Nn); %% fixed (min) level of sensor space variance


%==========================================================================
% Step 1: Optimise spatial priors over subjects
%==========================================================================

% create source components (Qp)
%==========================================================================
allind=[];
switch(type)
    
    case {'MSP','GS','ARD'}
        % create MSP spatial basis set in source space
        %------------------------------------------------------------------
        Qp    = {};
        LQpL  = {};
        %Ip    = ceil((1:Np)*Ns/Np);        % "random" selection of patches
        for i = 1:Np
            % Patch locations determined by Ip
            %--------------------------------------------------------------
            q               = QG(:,Ip(i));
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = UL*q;     
        end
        
%     case {'EBB'}
%         % create beamforming prior. See:
%         % Source reconstruction accuracy of MEG and EEG Bayesian inversion approaches.
%         %Belardinelli P, Ortiz E, Barnes G, Noppeney U, Preissl H. PLoS One. 2012;7(12):e51985.
%         %------------------------------------------------------------------
%         InvCov = spm_inv(YY);
%         allsource = zeros(Ns,1);
%         Sourcepower = zeros(Ns,1);
%         for bk = 1:Ns
%             normpower = 1/(UL(:,bk)'*UL(:,bk));
%             Sourcepower(bk) = 1/(UL(:,bk)'*InvCov*UL(:,bk));
%             allsource(bk) = Sourcepower(bk)./normpower;
%         end
%         allsource = allsource/max(allsource);   % Normalise
%         
%         Qp{1} = diag(allsource);
%         LQpL{1} = UL*diag(allsource)*UL';
   case {'EBB'}
        % create SMOOTH beamforming prior. 
        disp('NB smooth EBB algorithm !');
        %------------------------------------------------------------------
        InvCov = spm_inv(AYYA);
        allsource = sparse(Ns,1);
        Sourcepower = sparse(Ns,1);
        for bk = 1:Ns
            q               = QG(:,bk);
            
            smthlead = UL*q;     %% THIS IS WHERE THE SMOOTHNESS GETS ADDED
            normpower = 1/(smthlead'*smthlead);
            Sourcepower(bk) = 1/(smthlead'*InvCov*smthlead);
            allsource(bk) = Sourcepower(bk)./normpower;
        end
        allsource = allsource/max(allsource);   % Normalise
        
        Qp{1} = diag(allsource);
        LQpL{1} = UL*diag(allsource)*UL';
        
    case {'EBBgs'}  % NEW BEAMFORMER PRIOR!!
        % create beamforming prior- Juan David- Martinez Vargas
        %------------------------------------------------------------------
        allsource = zeros(Ntrials,Ns);
        for ii = 1:Ntrials
            InvCov = spm_inv(YYep{ii});
            Sourcepower = zeros(Ns,1);
            for bk = 1:Ns
                normpower = 1/(UL(:,bk)'*UL(:,bk));
                Sourcepower(bk) = 1/(UL(:,bk)'*InvCov*UL(:,bk));
                allsource(ii,bk) = Sourcepower(bk)./normpower;
            end
            
            Qp{ii}.q = allsource(ii,:);
        end
             
    case {'LOR','COH'}
        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = UL*UL';
        
        % add smoothness component in source space
        %------------------------------------------------------------------
        Qp{2}   = QG;
        LQpL{2} = UL*Qp{2}*UL';
        
    case {'IID','MMN'}
        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = UL*UL';
end

fprintf('Using %d spatial source priors provided\n',length(Qp));




% Inverse solution
%==========================================================================
QP     = {};
LQP    = {};
LQPL   = {};

% Get source-level priors (using all subjects)
%--------------------------------------------------------------------------
switch(type)
    
    case {'MSP','GS','EBBgs'}
        % Greedy search over MSPs
        %------------------------------------------------------------------
        Np    = length(Qp);
        Q     = zeros(Ns,Np); %% NB SETTING UP A NEW Q HERE
        for i = 1:Np
            Q(:,i) = Qp{i}.q;
        end
        Q = sparse(Q);
        
        % Multivariate Bayes (Here is performed the inversion)
        %------------------------------------------------------------------
        
        MVB   = spm_mvb(AY,UL,[],Q,Qe,16); %% Qe is identity with unit trace 
        
        % Accumulate empirical priors (New set of patches for the second inversion)
        %------------------------------------------------------------------
        % MVB.cp provides the final weights of the hyperparameters
        Qcp           = Q*MVB.cp;
        QP{end + 1}   = sum(Qcp.*Q,2);
        LQP{end + 1}  = (UL*Qcp)*Q';
        LQPL{end + 1} = LQP{end}*UL';
       
end

switch(type)
    
    case {'MSP','ARD'}
        
        % ReML / ARD   inversion
        %------------------------------------------------------------------
        
        
        %[Cy,h,Ph,F] = spm_sp_reml(AYYA,[],[Qe LQpL],1);
        [Cy,h,Ph,F] = spm_sp_reml(AYYA,[],[Qe LQpL],Nn);
        
        
        % Spatial priors (QP)
        %------------------------------------------------------------------
        % h provides the final weights of the hyperparameters
        Ne    = length(Qe);
        Np    = length(Qp);
        
        hp    = h(Ne + (1:Np));
        
        qp    = sparse(0);
        for i = 1:Np
            if hp(i) > max(hp)/128;
                qp  = qp + hp(i)*Qp{i}.q*Qp{i}.q';
            end
        end
        
        % Accumulate empirical priors (New set of patches for the second inversion)
        %------------------------------------------------------------------
        QP{end + 1}   = diag(qp);
        LQP{end + 1}  = UL*qp;
        LQPL{end + 1} = LQP{end}*UL';
end

switch(type)
    
    case {'IID','MMN','LOR','COH','EBB'}
        
        % or ReML - ARD (Here is performed the inversion)
        %------------------------------------------------------------------
        
        
        
        [Cy,h,Ph,F] = spm_reml_sc(AYYA,[],[Qe LQpL],Nn,-4,16,Q0);
        
        % Spatial priors (QP)
        %------------------------------------------------------------------
        % h provides the final weights of the hyperparameters
        Ne    = length(Qe);
        Np    = length(Qp);
        
        hp    = h(Ne + (1:Np));
        qp    = sparse(0);
        for i = 1:Np
            qp = qp + hp(i)*Qp{i};
        end
        
        % Accumulate empirical priors (New set of patches for the second inversion)
        %------------------------------------------------------------------
        QP{end + 1}   = diag(qp);
        LQP{end + 1}  = UL*qp;
        LQPL{end + 1} = LQP{end}*UL';
end


%==========================================================================
% Step 2: Re-estimate for each subject separately (fusing all modalities)
%==========================================================================

fprintf('Inverting subject 1\n')



% re-do ReML (with informative hyperpriors)
%----------------------------------------------------------------------
Np    = length(LQPL);       % Final number of priors
Ne    = length(Qe);         % Sensor noise prior


Q     = [{Q0} LQPL]; %% sensor corvariance prior:  Qe is identity with unit trace, LQPL is in the units of data

if rank(AYYA)~=size(A,1),
    rank(AYYA)
    size(AYYA,1)
    warning('AYYA IS RANK DEFICIENT');
end;


[Cy,h,Ph,F]= spm_reml_sc(AYYA,[],Q,Nn,-4,16,Q0);



%% recalculate F here

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


% accuracy; signal to noise (over sources)
%======================================================================
R2   = 100*(SST - SSR)/SST;
fprintf('Percent variance explained %.2f (%.2f)\n',full(R2),full(R2*VE));

% Save results
% DEMO: WARNING! These results are not coincident in format with
%                those generated in the SPM8
%======================================================================
inverse.type   = type;                 % inverse model
inverse.smooth = s;                    % smoothing coefficient
inverse.M      = M;                    % MAP projector (reduced)
inverse.J   = J;                    % Conditional expectation
inverse.Y      = Y;                    % ERP data (reduced)
inverse.L      = UL;                   % Lead-field (reduced)
inverse.qC     = Cq;                   % spatial covariance
inverse.tempU  = U;                    % temporal SVD 
inverse.V      = V;                    % temporal modes
inverse.qV     = Vq;                   % temporal correlations
inverse.T      = S;                    % temporal projector
inverse.U      = {A};                    % spatial projector
inverse.Is     = Is;                   % Indices of active dipoles
inverse.It     = It;                   % Indices of time bins
inverse.Ik     =Ik;                    %% indices of trials used
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
inverse.Ip=Ip;                          %% patch locations

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
