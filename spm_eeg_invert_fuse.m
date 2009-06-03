function [D] = spm_eeg_invert_fuse(varargin)
% ReML fusion of multiple forward models for EEG-MEG
% FORMAT [D] = spm_eeg_invert_fuse(D, val)
% ReML estimation of regularisation hyperparameters using the
% spatiotemporal hierarchy implicit in EEG data
% Requires:
% D{i}.inv{val}.inverse:
%
%     inverse.trials - D.events.types to invert
%     inverse.smooth - smoothness of source priors (0 to 1)
%     inverse.Np     - number of sparse priors per hemisphere
%     inverse.Nm     - maximum number of channel modes
%     inverse.type   - 'GS'  Greedy search on MSPs
%                      'ARD' ARD search on MSPs
%                      'MSP' GS and ARD multiple sparse priors
%                      'LOR' LORETA-like model
%                      'IID' LORETA and minimum norm
%     inverse.xyz    - (n x 3) locations of spherical VOIs
%     inverse.rad    - radius (mm) of VOIs
%     inverse.lpf    - band-pass filter - low -frequency cut-off (Hz)
%     inverse.hpf    - band-pass filter – high-frequency cut-off (Hz)
%     inverse.Lap    - switch for Laplace transform
%     inverse.sdv    - standard deviations of Gaussian temporal correlation
%     inverse.Han    - switch for Hanning window
%     inverse.Na     - number of most energetic dipoles
%     inverse.woi    - time window of interest ([start stop] in ms)
%
% Evaluates:
%
%     inverse.M      - MAP projector (reduced)
%     inverse.J      - Conditional expectation
%     inverse.L      - Lead field (reduced)
%     inverse.R      - Re-referencing matrix
%     inverse.qC     - spatial  covariance
%     inverse.qV     - temporal correlations
%     inverse.T      - temporal subspace
%     inverse.U      - spatial  subspace
%     inverse.Is     - Indices of active dipoles
%     inverse.Nd     - number of dipoles
%     inverse.pst    - peristimulus time
%     inverse.dct    - frequency range
%     inverse.F      - log-evidence
%     inverse.R2     - variance accounted for (%)
%     inverse.VE     - variance explained by temporal subspace
%     inverse.h      - covariance hyperparameters
%
%__________________________________________________________________________
% 
% This routine solves ill-posed linear models of the following form
%
%             y{1,...,t}  = L(1} * J{1,...,t}   +  e{{1,...,t}}
%             y{2,...,t}  = L(2}                   e{{2,...,t}}
%                  .
%                  .
%                  .
%             y{n,...,t}  = L(n}                   e{{n,...,t}}
%
% Under empirical priors on J{1,...,t} for n data-sets with t trial types.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_eeg_invert_fuse.m 3177 2009-06-03 08:47:41Z vladimir $
                        
 
% D - SPM MEEG object
%==========================================================================
[D, val]    = spm_eeg_inv_check(varargin{:}); 
[mod, list] = modality(D, 1, 1);
Nl          = numel(list);                       % number of forward models
inverse     = D.inv{val}.inverse;
 
% defaults
%--------------------------------------------------------------------------
try, type  = inverse.type;   catch, type  = 'GS';              end
try, s     = inverse.smooth; catch, s     = 0.6;               end
try, Np    = inverse.Np;     catch, Np    = 256;               end
try, Nm    = inverse.Nm;     catch, Nm    = [];                end
try, Nr    = inverse.Nr;     catch, Nr    = [];                end
try, NmTOL = inverse.NmTOL;  catch, NmTOL = -16;               end
try, NrTOL = inverse.NrTOL;  catch, NrTOL = 1;                 end
try, xyz   = inverse.xyz;    catch, xyz   = [0 0 0];           end
try, rad   = inverse.rad;    catch, rad   = 128;               end
try, lpf   = inverse.lpf;    catch, lpf   = 1;                 end
try, hpf   = inverse.hpf;    catch, hpf   = 256;               end
try, sdv   = inverse.sdv;    catch, sdv   = 4;                 end
try, Han   = inverse.Han;    catch, Han   = 1;                 end
try, Na    = inverse.Na;     catch, Na    = 1024;              end
try, woi   = inverse.woi;    catch, woi   = [];                end
 
try
    trial = inverse.trials;
catch
    trial = D.condlist;
end

 
%==========================================================================
% Spatial parameters
%==========================================================================
 
% Load Gain or Lead-field matrices
%--------------------------------------------------------------------------

for i = 1:Nl
    Ic{i}     = setdiff(meegchannels(D, list{i}), D.badchannels);
    [L{i}, D] = spm_eeg_lgainmat(D, [], D.chanlabels(Ic{i}));
    Nc(i)     = size(L{i},1);                      % number of channels
    Nd(i)     = size(L{i},2);                      % number of dipoles
end
 
% check the source space is the same
%--------------------------------------------------------------------------
if any(diff(Nd))
    warndlg('please ensure the number of dipoles is the same')
    return
end
 
 
% assume radii are the same for all VOI
%--------------------------------------------------------------------------
Nd    = Nd(1);                                     % number of dipoles
Nv    = size(xyz,1);                               % number of VOI
if length(rad) ~= Nv
    rad = rad(1)*ones(Nv,1);
else
    rad = rad(:);
end
 
% Restrict source space
%--------------------------------------------------------------------------
vert  = D.inv{val}.mesh.tess_mni.vert;
face  = D.inv{val}.mesh.tess_mni.face;
Is    = sparse(Nd,1);
for i = 1:Nv
    Iv = sum([vert(:,1) - xyz(i,1), ...
              vert(:,2) - xyz(i,2), ...
              vert(:,3) - xyz(i,3)].^2,2) < rad(i)^2;
    Is = Is | Iv;
end
Is    = find(Is);
vert  = vert(Is,:);
for i = 1:Nl
    L{i} = L{i}(:,Is);
end
Ns    = length(Is);
 
 
% Compute spatial coherence: Diffusion on a normalised graph Laplacian GL
%==========================================================================
 
fprintf('Computing Green''s function from graph Laplacian:')
%--------------------------------------------------------------------------
A     = spm_eeg_inv_meshdist(vert,face,0);
GL    = A - spdiags(sum(A,2),0,Nd,Nd);
GL    = GL*s/2;
Qi    = speye(Nd,Nd);
QG    = sparse(Nd,Nd);
for i = 1:8
    QG = QG + Qi;
    Qi = Qi*GL/i;
end
clear Qi
QG    = QG.*(QG > exp(-8));
QG    = QG*QG;
QG    = QG(Is,Is);
fprintf(' - done\n')
 
 
% Project to channel modes (U)
%--------------------------------------------------------------------------
if isempty(Nm); 
  Nm = zeros(1,Nl); 
elseif length(Nm) == 1; 
  Nm = ones(1,Nl)*Nm; 
end
for i = 1:Nl
    if Nm(i) > Nc(i)
      error(sprintf('Specified Nm (%d) must be less than Num Channels (%d)',Nm(i),Nc(i)));
    end
    if ~Nm(i)
        U{i}  = spm_svd(L{i}*L{i}',exp(NmTOL));
        Nm(i) = size(U{i},2);
    else
        U{i}  = spm_svd(L{i}*L{i}',exp(-32));
        U{i}  = U{i}(:,1:Nm(i));
    end
    UL{i} = U{i}'*L{i};
    fprintf('Using %i spatial modes\n',Nm(i))
end

%==========================================================================
% Temporal parameters
%==========================================================================
Nt     = length(trial);                % number of trial types
if isempty(Nr);
    Nr = zeros(1,Nl);
elseif length(Nr) ~= Nl
    Nr = ones(1,Nl)*Nr(1);
end

% Peristimulus time-window 
%--------------------------------------------------------------------------
w = [D.time(1) D.time(end)]*1000;

if ~isempty(woi)
    woi(1) = max(woi(1),w(1));
    woi(2) = min(woi(2),w(2));
else
    woi    = w;
end
 
% get time bins to sample
%--------------------------------------------------------------------------
It    = D.indsample(woi/1000);
Nb    = fix(It(end) - It(1) + 1);              % number of time bins
It    = fix(linspace(It(1),It(2),Nb));


% Peri-stimulus time
%--------------------------------------------------------------------------
woi   = fix(woi);
pst   = 1000*D.time(It);                       % peristimulus time (ms)
dur   = (pst(end) - pst(1))/1000;              % duration (s)
dct   = ([1:Nb] - 1)/2/dur;                    % DCT frequencies (Hz)
 
% Serial correlations
%----------------------------------------------------------------------
K     = exp(-(pst - pst(1)).^2/(2*sdv^2));
K     = toeplitz(K);
qV    = sparse(K*K');
 
%  Frequency subspace
%--------------------------------------------------------------------------
T     = spm_dctmtx(Nb,Nb);
j     = find( (dct >= lpf) & (dct <= hpf) );
T     = T(:,j);
dct   = dct(j);

% get data (with temporal filtering)
%==========================================================================
UY    = {};
AY    = {};
for i = 1:Nl
    
    if Nr(i) > size(T,2)
       warning(sprintf(['Specified Nr (%d) greater than Num Freq components' ...
            ' between %d and %d Hz; so using %d'],Nr,lpf,hpf,size(T,2)));
       Nr(i) = size(T,2);
    end
    
    % get channels
    %----------------------------------------------------------------------
    for j = 1:Nt
        UY{i,j} = sparse(0);
        c     = D.pickconditions(trial{j});
        Ne    = length(c);
        for k = 1:Ne
            UY{i,j} = UY{i,j} + U{i}'*squeeze(D(Ic{i},It,c(k)))*T/Ne;
        end
    end
end

% rescale under i.i.d. priors
%--------------------------------------------------------------------------
UL1   = UL{1}'*inv(UL{1}*UL{1}');
for i = 1:Nl
    for j = 1:Nt
        scale(i,j) = sign(trace(UY{i,j}'*UL{i}*UL1*UY{1,j}));
        scale(i,j) = scale(i,j)/sqrt(trace(UY{i,j}'*UY{i,j})/Nm(i));
    end
end
scale = mean(scale,2);
for i = 1:Nl
    for j = 1:Nt
        UY{i,j} = UY{i,j}*scale(i);
    end
    UL{i} = UL{i}/sqrt(trace(UL{i}*UL{i}')/Nm(i));
end
 
% and apply Hanning operator (if requested) 
%--------------------------------------------------------------------------
if Han
    W = T'*sparse(1:Nb,1:Nb,spm_hanning(Nb))*T;
else
    W = 1;
end
for i = 1:Nl
    WY{i} = UY(i,:);
    WY{i} = spm_cat(WY{i}(:))*W;
end
 
% temporal projector (at most 8 modes) S = T*V
%==========================================================================
for i = 1:Nl
    if ~Nr(i)
        [v u] = spm_svd(WY{i}',NrTOL);            % temporal modes
        Nr(i) = size(v,2);
    else
        [v u] = spm_svd(WY{i}',exp(-16));
        v     = v(:,       1:Nr(i));
        u     = u(1:Nr(i), 1:Nr(i));
    end
    VE(i) = sum(sum(u.^2))/sum(sum(WY{i}.^2));
    V{i}  = v; 
    fprintf('Using %i temporal modes, accounting for %0.2f percent variance\n',Nr(i),full(100*VE(i)))
end

V     = orth(full(spm_cat(V)));
Nr    = sum(Nr);
VE    = mean(VE); 
fprintf('Using %i temporal modes in total...\n',Nr)
fprintf('...accounting for %0.2f percent variance on average\n',full(100*VE))

% projection and whitening (noting that UY = U*Y*T)
%--------------------------------------------------------------------------
S     = T*V;                                     % temporal projector
iV    = inv(S'*qV*S);                            % precision (mode)
Vq    = S*iV*S';                                 % precision (time)
for i = 1:Nl
    for j = 1:Nt
        UY{i,j} = UY{i,j}*V;
        AY{i,j} = UY{i,j}*sqrtm(iV);
    end
end
 
% concatenate data
%--------------------------------------------------------------------------
UY    = spm_cat(UY,1);
AY    = spm_cat(AY);
YY    = AY*AY';
G     = spm_cat(UL(:));
 
% create sensor components (Qe); one for each modality
%==========================================================================
for i = 1:Nl
    QE{i,i} = sparse(Nm(i),Nm(i));
end
for i = 1:Nl
    Q      = QE;
    Q{i,i} = U{i}'*U{i};
    Qe{i}  = spm_cat(Q);
end
Ne    = length(Qe);
 
 
% create source components (Qp)
%==========================================================================
switch(type)
 
    case {'MSP','GS','ARD'}
 
        % create MSP spatial basis set in source space
        %------------------------------------------------------------------
        Qp    = {};
        LQpL  = {};
        Ip    = ceil([1:Np]*Ns/Np);
        for i = 1:Np
 
            % left hemisphere
            %--------------------------------------------------------------
            q               = QG(:,Ip(i));
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = G*q;
 
            % right hemisphere
            %--------------------------------------------------------------
            [d j] = min(sum([vert(:,1) + vert(Ip(i),1), ...
                vert(:,2) - vert(Ip(i),2), ...
                vert(:,3) - vert(Ip(i),3)].^2,2));
            q               = QG(:,j);
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = G*q;
 
            % bilateral
            %--------------------------------------------------------------
            q               = QG(:,Ip(i)) + QG(:,j);
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = G*q;
 
        end
 
    case {'LOR','COH'}
 
        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = G*G';
 
        % add smoothness component in source space
        %------------------------------------------------------------------
        Qp{2}   = QG;
        LQpL{2} = G*Qp{2}*G';
 
 
    case {'IID','MMN'}
 
        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = G*G';
 
end
 
% Inverse solution
%==========================================================================
QP     = {};
LQP    = {};
LQPL   = {};
 
% get source-level priors (using all subjects)
%--------------------------------------------------------------------------
switch(type)
 
    case {'MSP','GS'}
 
        % Greedy search over MSPs
        %------------------------------------------------------------------
        Np    = length(Qp);
        Q     = sparse(Ns,Np);
        for i = 1:Np
            Q(:,i) = Qp{i}.q;
        end
 
        % Multivariate Bayes
        %------------------------------------------------------------------
        MVB   = spm_mvb(AY,G,[],Q,Qe,16);
        F1    = MVB.F;
 
        % Accumulate empirical priors
        %------------------------------------------------------------------
        Qcp           = Q*MVB.cp;
        QP{end + 1}   = sum(Qcp.*Q,2);
        LQP{end + 1}  = (G*Qcp)*Q';
        LQPL{end + 1} = LQP{end}*G';

end
 
switch(type)
 
    case {'MSP','ARD'}
 
        % or ReML - ARD
        %------------------------------------------------------------------
        Q            = {Qe{:} LQpL{:}};
        [Cy,h,Ph,F1] = spm_sp_reml(YY,[],Q,Nr*Nt);

        % Spatial priors (QP)
        %------------------------------------------------------------------
        Ne    = length(Qe);
        Np    = length(Qp);
        hp    = h([1:Np] + Ne);
        qp    = sparse(0);
        for i = 1:Np
            if hp(i) > max(hp)/128;
                qp = qp + hp(i)*Qp{i}.q*Qp{i}.q';
            end
        end
 
        % Accumulate empirical priors
        %------------------------------------------------------------------
        QP{end + 1}   = diag(qp);
        LQP{end + 1}  = G*qp;
        LQPL{end + 1} = LQP{end}*G';
 
end

switch(type)
 
    case {'IID','MMN','LOR','COH'}
 
        % or ReML - ARD
        %------------------------------------------------------------------
        Q     = {Qe{:} LQpL{:}};
        [Cy,h,Ph,F1] = spm_reml_sc(YY,[],Q,Nr*Nt);

        % Spatial priors (QP)
        %------------------------------------------------------------------
        Ne    = length(Qe);
        Np    = length(Qp);
        hp    = h([1:Np] + Ne);
        qp    = sparse(0);
        for i = 1:Np
            qp = qp + hp(i)*Qp{i};
        end
 
        % Accumulate empirical priors
        %------------------------------------------------------------------
        QP{end + 1}   = diag(qp);
        LQP{end + 1}  = G*qp;
        LQPL{end + 1} = LQP{end}*G';
 
end
 
% re-estimate
%==========================================================================
 
% using spatial priors from group analysis
%--------------------------------------------------------------------------
Np    = length(LQPL);
Q     = {Qe{:} LQPL{:}};
 
% re-do ReML
%--------------------------------------------------------------------------
[Cy,h,Ph,F] = spm_reml_sc(YY,[],Q,Nr*Nt);

% Data ID
%==========================================================================
ID  = spm_data_id(YY);
 
% Covariance: sensor space - Ce and source space - L*Cp
%--------------------------------------------------------------------------
Cp    = sparse(0);
LCp   = sparse(0);
hp    = h([1:Np] + Ne);
for j = 1:Np
     Cp =  Cp + hp(j)*QP{j};
    LCp = LCp + hp(j)*LQP{j};
end

 
% MAP estimates of instantaneous sources
%==========================================================================
iC    = inv(Cy);
M     = LCp'*iC;
 
% conditional covariance (leading diagonal)
% Cq    = Cp - Cp*L'*iC*L*Cp;
%--------------------------------------------------------------------------
Cq    = Cp - sum(LCp.*M')';
 
% evaluate conditional expectation (of the sum over trials)
%--------------------------------------------------------------------------
for j = 1:Nt
 
    % trial-type specific source reconstruction
    %----------------------------------------------------------------------
    J{j}  = M*UY{j};
 
    % sum of squares
    %----------------------------------------------------------------------
    SSR(j)   = sum(var((UY{j} - G*J{j}),0,2));
    SST(j)   = sum(var(UY{j},0,2));
 
end
 
% assess accuracy; signal to noise (over sources)
%==========================================================================
R2   = 100*(SST - SSR)./SST
R2   = mean(R2);
fprintf('Percent variance explained %.2f (%.2f)\n',full(R2),full(R2*VE))
 
% Save results
%==========================================================================
inverse.type   = type;                 % inverse model
inverse.smooth = s;                    % smoothness (0 - 1)
inverse.xyz    = xyz;                  % VOI (XYZ)
inverse.rad    = rad;                  % VOI (radius)
inverse.scale  = scale;                % scale factors
inverse.M      = M;                    % MAP projector (reduced)
inverse.J      = J;                    % Conditional expectation
inverse.Y      = UY;                   % ERP data (reduced)
inverse.L      = G;                    % Lead-field (reduced)
inverse.R      = 1;                    % Re-referencing matrix
inverse.qC     = Cq;                   % spatial covariance
inverse.qV     = Vq;                   % temporal correlations
inverse.T      = S;                    % temporal subspace
inverse.U      = spm_cat(diag(U));     % spatial subspace
inverse.Is     = Is;                   % Indices of active dipoles
inverse.It     = It;                   % Indices of time bins
inverse.Ic     = Ic;                   % Indices of good channels
inverse.Nd     = Nd;                   % number of dipoles
inverse.pst    = pst;                  % peristimulus time
inverse.dct    = dct;                  % frequency range
inverse.F1     = F1;                   % log-evidence (first step; with MSP complexity)
inverse.F      = F;                    % log-evidence (second step; MSP complexity ignored)
inverse.ID     = ID;                   % data ID
inverse.R2     = R2;                   % variance accounted for (%)
inverse.VE     = VE;                   % variance explained
inverse.woi    = woi;                  % time-window inverted
inverse.h      = h;                    % covariance hyperparameters (second step)
inverse.Nm     = Nm;                   % number of spatial modes
inverse.Nr     = Nr;                   % number of spatial modes

%--------------------------------------------------------------------------

D.inv{val}.inverse = inverse;
D.inv{val}.method  = 'Imaging';

% and delete old contrasts
%--------------------------------------------------------------------------
try
    D.inv{val} = rmfield(D.inv{val},'contrast');
end
 
% display
%==========================================================================
spm_eeg_invert_display(D);
drawnow

