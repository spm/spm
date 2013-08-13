function [D] = spm_eeg_invert_classic(D,val);
%
%
% A trimmed down version of the spm_eeg_invert() routine
% see also spm_eeg_invert_noscale.m
%
% Created by:	Jose David Lopez - ralph82co@gmail.com
%				Gareth Barnes - g.barnes@fil.ion.ucl.ac.uk
%				Vladimir Litvak - litvak.vladimir@gmail.com
%
% Date: January 2012
%
% This version is for single subject single modality analysis and therefore
% contains none of the associated scaling factors.
% No symmetric priors are used in this implementation (just single patches)
% There is an option for a Beamforming prior : inversion type 'EBB'
%% Ip: determines the mesh vertices where the patches will be centred
%% Nm number of spatial modes to use
%% Nt number of temporal modes to use
%% SmthInit, default 0.6, spatial smoothing used to set size of patches
% Gareth
%% if SHUFFLEADS is no zero then this will randomly assign lead fields to different channels
%%
%%The code was used in
%% L�pez, J. D., Penny, W. D., Espinosa, J. J., Barnes, G. R. (2012).
% A general Bayesian treatment for MEG source reconstruction incorporating lead field uncertainty.
% Neuroimage 60(2), 1194-1204 doi:10.1016/j.neuroimage.2012.01.077.

% $Id: spm_eeg_invert_classic.m 5609 2013-08-13 16:48:17Z gareth $



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
try, hpf  = inverse.hpf;    catch, hpf  = 48;        end %% need to one day put these the correct way round
try, lpf  = inverse.lpf;    catch, lpf  = 0;       end
try, sdv  = inverse.sdv;    catch, sdv  = 4;        end
try, Han  = inverse.Han;    catch, Han  = 1;        end
try, woi  = inverse.woi;    catch, woi  = [];       end
try, Nm=inverse.Nm; catch, Nm=[];end;
try, Nt=inverse.Nt; catch, Nt=[];end; %% fixed number of temporal modes
try, Ip=inverse.Ip; catch Ip=[]; end;%
try SHUFFLELEADS=inverse.SHUFFLELEADS;catch SHUFFLELEADS=[];end;



% defaults
%--------------------------------------------------------------------------
type = inverse.type;	% Type of inversion scheme


% get specified modalities to invert (default to all)

%--------------------------------------------------------------------------
modalities = D.inv{val}.forward.modality;		% MEG in this case
Nmax  = 16;			% max number of temporal modes

% check lead fields and get number of dipoles (Nd) and channels (Nc)
%==========================================================================

fprintf('Checking leadfields')
[L D] = spm_eeg_lgainmat(D);	% Generate/load lead field
Nd=size(L,2);
if ~isempty(Ip)
    Np   = length(Ip);				% Number of priors/3 for GS, ARD, MSP
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
Nd    = size(L,2);		% Number of dipoles

fprintf(' - done\n')


% Compute spatial coherence: Diffusion on a normalised graph Laplacian GL
%==========================================================================

fprintf('Computing Green function from graph Laplacian:')
%--------------------------------------------------------------------------
vert  = D.inv{val}.mesh.tess_mni.vert;
face  = D.inv{val}.mesh.tess_mni.face;
A     = spm_mesh_distmtx(struct('vertices',vert,'faces',face),0);
GL    = A - spdiags(sum(A,2),0,Nd,Nd);
GL    = GL*s/2;				% Smoother
Qi    = speye(Nd,Nd);
QG    = sparse(Nd,Nd);
for i = 1:8					% Taylor series approximation
    QG = QG + Qi;
    Qi = Qi*GL/i;
end
QG    = QG.*(QG > exp(-8));		% Eliminate small values
QG    = QG*QG;				% Guarantee positive semidefinite matrix
%%QG=QG./sum(sum(QG)); %% normalise area
clear Qi A GL
fprintf(' - done\n')


% check for (e.g., empty-room) sensor components (in Qe)
%==========================================================================
QE = 1;						% No empty room


%==========================================================================
% Spatial projectors (adjusting for different Lead-fields)
%==========================================================================

fprintf('Optimising and aligning spatial modes ...\n')

% eliminate low SNR spatial modes
%------------------------------------------------------------------

if isempty(Nm),
    [U,ss,vv]    = spm_svd((L*L'),exp(-16));
    A     = U';					% spatial projector A
    UL    = A*L;
    
else
    [U,ss,vv]    = spm_svd((L*L'),0);
    if length(ss)<Nm,
        error('Not this many spatial modes in lead fields');
        disp('number available');
        length(ss)
    end;
        
    ss=ss(1:Nm);
    disp('using preselected number spatial modes !');
    A     = U(:,1:Nm)';					% spatial projector A
    UL    = A*L;
end;

Nm    = size(UL,1);			% Number of spatial projectors
% Plot spatial projectors
%------------------------------------------------------------------
% figure
% loglog(ss);					% Plot of singular values
% title('Spatial projector');
% xlabel('Eigenvalues');
% ylabel('Amplitude');
clear ss vv

% Report
%----------------------------------------------------------------------
fprintf('Using %d spatial modes',Nm)

% None dipole is eliminated
%--------------------------------------------------------------------------
Is    = 1:Nd;				% Accepted dipoles
Ns    = length(Is);			% Ns = Nd in this case


%==========================================================================
% Temporal projector
%==========================================================================

% Time-window of interest
%----------------------------------------------------------------------

if isempty(woi)
    w      = 1000*[min(D.time) max(D.time)];
else
    w=woi; %% in milliseconds
end;

It     = (w/1000 - D.timeonset)*D.fsample + 1;
It     = max(1,It(1)):min(It(end), length(D.time));
It	   = fix(It);

% Peristimulus time
%----------------------------------------------------------------------
pst    = 1000*D.time;					% peristimulus time (ms)
pst    = pst(It);						% windowed time (ms)
dur    = (pst(end) - pst(1))/1000;		% duration (s)
dct    = (It - It(1))/2/dur;			% DCT frequencies (Hz)
Nb     = length(It);					% number of time bins

% Serial correlations
%----------------------------------------------------------------------
K      = exp(-(pst - pst(1)).^2/(2*sdv^2));
K      = toeplitz(K);
qV     = sparse(K*K');

% Confounds and temporal subspace
%----------------------------------------------------------------------

T      = spm_dctmtx(Nb,Nb);			% use plot(T) here!

j      = find( (dct >= lpf) & (dct <= hpf) ); %% THis is the wrong way round but leave for nowfor compatibility with spm_eeg_invert
T      = T(:,j);					% Apply the filter to discrete cosines
dct    = dct(j);					% Frequencies accepted

%% Hanning window
%----------------------------------------------------------------------

if Han
    W  = sparse(1:Nb,1:Nb,spm_hanning(Nb)); %% use hanning unless specified
else
    W=1;
end;
% W  = 1;					% Apply Hanning if desired!

% get temporal covariance (Y'*Y) to find temporal modes
%======================================================================
% Note: The variable YY was replaced with YTY because it is
% duplicated in the original script, causing confusion.

Y      = A*D(Ic,It,1);	% Data samples in spatial modes (virtual sensors)
YTY    = Y'*Y;			% Covariance in temporal domain

% Apply any Hanning and filtering
%------------------------------------------------------------------
YTY         = W'*YTY*W;		% Hanning
YTY         = T'*YTY*T;		% Filter

% Plot temporal projectors
%------------------------------------------------------------------
% figure
% imagesc(YTY);		% Plot of frequency map
% title('Temporal projector');
% xlabel('Frequency (Hz)');
% ylabel('Frequency (Hz)');

% temporal projector (at most Nrmax modes) S = T*V
%======================================================================

if isempty(Nt),
    
    [U E]  = spm_svd(YTY,exp(-8));			% get temporal modes
    if isempty(U),
        warning('nothing found using spm svd, using svd');
        [U E]  = svd(YTY);			% get temporal modes
    end;
    E      = diag(E)/trace(YTY);			% normalise variance
    Nr     = min(length(E),Nmax);			% number of temporal modes
    Nr=max(Nr,1); %% use at least one mode
else
    [U E]  = svd(YTY);			% get temporal modes
    
    E      = diag(E)/trace(YTY);			% normalise variance
    disp('Fixed number of temporal modes');
    Nr=Nt;
end;
V      = U(:,1:Nr);						% temporal modes
VE     = sum(E(1:Nr));					% variance explained

fprintf('Using %i temporal modes, ',Nr)
fprintf('accounting for %0.2f percent average variance\n',full(100*VE))

% projection and whitening
%----------------------------------------------------------------------
S      = T*V;							% temporal projector
Vq     = S*pinv(S'*qV*S)*S';			% temporal precision


% get spatial covariance (Y*Y') for Gaussian process model
%======================================================================

% stack (scaled aligned data) over modalities
%--------------------------------------------------------------

Y    = A*D(Ic,It,1)*S;		% Load data in spatial and temporal modes

% accumulate first & second-order responses
%--------------------------------------------------------------
YY   = sparse(Y*Y');            % Data covariance in spatial mode

% generate sensor error components (Qe)
%==========================================================================

% assuming equal noise over subjects (Qe) and modalities AQ
%--------------------------------------------------------------------------
AQeA   = A*QE*A';			% Note that here it is A*A'
Qe{1}  = AQeA/(trace(AQeA)); % it means IID noise in virtual sensor space

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
        %Ip    = ceil((1:Np)*Ns/Np);		% "random" selection of patches
        for i = 1:Np
            % First set (Not left hemisphere)
            %--------------------------------------------------------------
            q               = QG(:,Ip(i));
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = UL*q;
            
            % 			% Extended set (add 256 priors)
            % 			%--------------------------------------------------------------
            % 			[dum,j] = min(sum([vert(:,1) + vert(Ip(i),1), ...
            % 				vert(:,2) - vert(Ip(i),2), ...
            % 				vert(:,3) - vert(Ip(i),3)].^2,2));
            % 			q               = QG(:,j);
            % 			Qp{end + 1}.q   = q;
            % 			LQpL{end + 1}.q = UL*q;
            %
            % 			% bilateral (add another 256 priors)
            % 			%--------------------------------------------------------------
            % 			% The bilateral patches are important with temporal
            % 			% lobe activity (synchronous sources)
            % 			q               = QG(:,Ip(i)) + QG(:,j);
            % 			Qp{end + 1}.q   = q;
            % 			LQpL{end + 1}.q = UL*q;
            
            
        end
        
    case {'EBB'}
        % create beamforming prior. See:
        % Source reconstruction accuracy of MEG and EEG Bayesian inversion approaches. 
        %Belardinelli P, Ortiz E, Barnes G, Noppeney U, Preissl H. PLoS One. 2012;7(12):e51985. 
        %------------------------------------------------------------------
        InvCov = spm_inv(YY);
        allsource = zeros(Ns,1);
        Sourcepower = zeros(Ns,1);
        for bk = 1:Ns
            normpower = 1/(UL(:,bk)'*UL(:,bk));
            Sourcepower(bk) = 1/(UL(:,bk)'*InvCov*UL(:,bk));
            allsource(bk) = Sourcepower(bk)./normpower;
        end
        allsource = allsource/max(allsource);	% Normalise
        
        Qp{1} = diag(allsource);
        LQpL{1} = UL*diag(allsource)*UL';
        
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
    
    case {'MSP','GS'}
        % Greedy search over MSPs
        %------------------------------------------------------------------
        Np    = length(Qp);
        Q     = zeros(Ns,Np);
        for i = 1:Np
            Q(:,i) = Qp{i}.q;
        end
        Q = sparse(Q);
        
        % Multivariate Bayes (Here is performed the inversion)
        %------------------------------------------------------------------
        MVB   = spm_mvb(Y,UL,[],Q,Qe,16);
        
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
        
        % ReML - ARD (Here is performed the inversion)
        %------------------------------------------------------------------
        [Cy,h,Ph,F] = spm_sp_reml(YY,[],[Qe LQpL],1);
        
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
        Q0          = exp(-2)*trace(YY)*Qe{1}/trace(Qe{1});
        [Cy,h,Ph,F] = spm_reml_sc(YY,[],[Qe LQpL],1,-4,16,Q0);
        
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

% generate sensor component (Qe) per modality
%----------------------------------------------------------------------
AQeA  = A*QE*A';				% Again it is A*A'
AQ    = AQeA/(trace(AQeA));


% using spatial priors from group analysis
%----------------------------------------------------------------------
Np    = length(LQPL);		% Final number of priors
Ne    = length(Qe);			% Sensor noise prior
Q     = [Qe LQPL];

% re-do ReML (with informative hyperpriors)
% Here is performed the second inversion
%======================================================================
Q0          = exp(-2)*trace(YY)*AQ/trace(AQ);
[Cy,h,Ph,F] = spm_reml_sc(YY,[],Q,1,-4,16,Q0);

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
J = M*Y;

% sum of squares
%------------------------------------------------------------------
SSR  = sum(var((Y - UL*J),0,2));
SST  = sum(var(Y,0,2));

% accuracy; signal to noise (over sources)
%======================================================================
R2   = 100*(SST - SSR)/SST;
fprintf('Percent variance explained %.2f (%.2f)\n',full(R2),full(R2*VE));

% Save results
% DEMO: WARNING! These results are not coincident in format with
%                those generated in the SPM8
%======================================================================
inverse.type   = type;                 % inverse model
inverse.smooth = s;                    % smoothness (0 - 1)
inverse.M      = M;                    % MAP projector (reduced)
inverse.J{1}   = J;                    % Conditional expectation
inverse.Y      = Y;                    % ERP data (reduced)
inverse.L      = UL;                   % Lead-field (reduced)
inverse.qC     = Cq;                   % spatial covariance
inverse.qV     = Vq;                   % temporal correlations
inverse.T      = S;                    % temporal projector
inverse.U      = {A};                    % spatial projector
inverse.Is     = Is;                   % Indices of active dipoles
inverse.It     = It;                   % Indices of time bins
inverse.Ic     = {Ic};                   % Indices of good channels
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
