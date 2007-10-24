function [D] = spm_eeg_invert_single(D)
% ReML inversion of forward model for EEG-EMG
% FORMAT [D] = spm_eeg_invert(D)
% ReML estimation of regularisation hyperparameters using the 
% spatio-temporal hierarchy implicit in EEG data
% Requires:
% D.inv{val}.inverse:
%
%     inverse.trials - D.events.types to invert
%     inverse.smooth - smoothness of source priors (0 to 1)
%     inverse.Np     - number of sparse priors per hemisphere
%     inverse.Nm     - maximum number of channel modes
%     inverse.type   - 'MSP' multiple sparse priors
%                      'LOR' LORETA-like model
%                      'IID' LORETA and minimum norm
%     inverse.xyz    - (n x 3) locations of spherical VOIs
%     inverse.rad    - radius (mm) of VOIs
%     inverse.lpf    - band-pass filter - low  frequency cutoff (Hz)
%     inverse.hpf    - band-pass filter - high frequency cutoff (Hz)
%     inverse.Lap    - switch for Laplace transform
%     inverse.sdv    - standard devations of Gaussian temporal correlation
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
%     inverse.pst    - pers-stimulus time
%     inverse.dct    - frequency range
%     inverse.F      - log-evidence
%     inverse.R2     - variance accounted for (%)
%__________________________________________________________________________
% Modified by Rik Henson to allow more arguments to be passed       4/6/07

% D - SPM data structure
%==========================================================================
model = D.inv{D.val};

% defaults
%--------------------------------------------------------------------------
try, trial = model.inverse.trials; catch, trial = D.events.types; end
try, type  = model.inverse.type;   catch, type  = 'GS';           end
try, s     = model.inverse.smooth; catch, s     = 0.6;            end
try, Np    = model.inverse.Np;     catch, Np    = 256;            end
try, Nm    = model.inverse.Nm;     catch, Nm    = 128;             end
try, xyz   = model.inverse.xyz;    catch, xyz   = [0 0 0];        end
try, rad   = model.inverse.rad;    catch, rad   = 128;            end
try, lpf   = model.inverse.lpf;    catch, lpf   = 1;              end
try, hpf   = model.inverse.hpf;    catch, hpf   = 256;            end
try, sdv   = model.inverse.sdv;    catch, sdv   = 4;              end
try, Han   = model.inverse.Han;    catch, Han   = 1;              end
try, Na    = model.inverse.Na;     catch, Na    = 1024;           end
try, woi   = model.inverse.woi;    catch, woi   = [];             end

% Load Gain or Lead field matrix
%--------------------------------------------------------------------------
try
    L     = load(model.forward.gainmat);
catch
    [p f] = fileparts(model.forward.gainmat);
    L     = load(f);
end
name  = fieldnames(L);
L     = sparse(getfield(L, name{1}));

% Time-window of interest
%--------------------------------------------------------------------------
if isempty(woi) 
    woi   = round([-D.events.start D.events.stop]*1000/D.Radc);
end
It    = round(woi*(D.Radc/1000)) + D.events.start;
It    = max(1,It(1)):min(It(end),size(D.data,2));

% parameters
%==========================================================================
Nb    = length(It);                              % number of time bins
Nc    = size(L,1);                               % number of channels
Nd    = size(L,2);                               % number of dipoles
Nv    = size(xyz,1);                             % number of VOI

% assume radii are the same for all VBOI
%--------------------------------------------------------------------------
if length(rad) ~= Nv
    rad = rad(1)*ones(Nv,1);
else
    rad = rad(:);
end

% Peri-stimulus time
%--------------------------------------------------------------------------
pst   = (It - D.events.start - 1)/D.Radc*1000;   % peristimulus time (ms)
dur   = (pst(end) - pst(1))/1000;                % duration (s)
dct   = (It - It(1))/2/dur;                      % DCT frequenices (Hz)

% Serial correlations
%--------------------------------------------------------------------------
K     = exp(-(pst - pst(1)).^2/(2*sdv^2));
K     = toeplitz(K);
qV    = sparse(K*K');

% Confounds and temporal subspace
%--------------------------------------------------------------------------
T     = spm_dctmtx(Nb,Nb);
i     = (dct > lpf) & (dct < hpf);
T     = T(:,i);
dct   = dct(i);

% get data (with temporal reduction)
%==========================================================================
Ic    = setdiff(D.channels.eeg, D.channels.Bad);
Nt    = length(trial);
for i = 1:Nt
    Y{i}  = sparse(0); 
    if isfield(D.events,'reject')
        c = find(D.events.code == trial(i) & ~D.events.reject);
    else
 	    c = find(D.events.code == trial(i));
    end
    for j = 1:length(c)
        Y{i} = Y{i} + squeeze(D.data(Ic,It,c(j)))*T;
    end
end

% temporal covariance (with Hanning if requested)
%--------------------------------------------------------------------------
YY    = sparse(0);
if Han
    W = T'*diag(spm_hanning(Nb))*T;
else
    W = T'*T;
end
for i = 1:Nt
    YY  = YY + W'*Y{i}'*Y{i}*W;
end

% eliminate unnecessary temporal modes
%--------------------------------------------------------------------------
S     = spm_svd(YY,exp(-4));
T     = T*S;
Nr    = size(T,2);                               % number of temporal modes
iV    = inv(T'*qV*T);                            % precision (mode space)
Vq    = T*iV*T';
fprintf('Using %i temporal modes\n',Nr)


% Project onto temporal modes (S)
%--------------------------------------------------------------------------
YY    = sparse(0);
for i = 1:Nt
    Y{i} = Y{i}*S;
    YY   = YY + Y{i}*iV*Y{i}';
end
UY    = spm_cat(Y)*kron(speye(Nt,Nt),sqrtm(iV));


% Project to channel modes (U)
%--------------------------------------------------------------------------
U     = spm_svd(L*L',exp(-32));
try
    U = U(:,1:Nm);
end
Nm    = size(U,2);
UY    = U'*UY;
YY    = U'*YY*U;
L     = U'*L;
fprintf('Using %i spatial modes\n',Nm)

% Restrict source space
%==========================================================================
vert  = model.mesh.tess_mni.vert;
face  = model.mesh.tess_mni.face;
Is    = sparse(Nd,1);
for i = 1:Nv
        Iv = sum([vert(:,1) - xyz(i,1), ...
                  vert(:,2) - xyz(i,2), ...
                  vert(:,3) - xyz(i,3)].^2,2) < rad(i)^2;
        Is = Is | Iv;
end
Is    = find(Is);
vert  = vert(Is,:);
L     = L(:,Is);
Ns    = length(Is);
Na    = min(Na,Ns);

% Compute spatial coherence: Diffusion on a normalised graph Laplacian GL
%==========================================================================

if ~strcmp(type,'IID')
    
    fprintf('Computing Green''s function from graph Laplacian:')
    %----------------------------------------------------------------------
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
    
end


% covariance components
%==========================================================================

% sensor noise (accommodating re-reference)
%--------------------------------------------------------------------------
Qe{1} = U'*U;

% create source compeonts
%--------------------------------------------------------------------------
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
            LQpL{end + 1}.q = L*q;
            
            % right hemisphere
            %--------------------------------------------------------------
            [d j] = min(sum([vert(:,1) + vert(Ip(i),1), ...
                             vert(:,2) - vert(Ip(i),2), ...
                             vert(:,3) - vert(Ip(i),3)].^2,2));
            q               = QG(:,j);
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = L*q;
            
            % bilateral
            %--------------------------------------------------------------
            q               = QG(:,Ip(i)) + QG(:,j);
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = L*q;
            
        end

    case {'LOR','COH'}

        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = L*L';

        % add smoothness component in source space
        %------------------------------------------------------------------
        Qp{2}   = QG;
        LQpL{2} = L*Qp{2}*L';


    case {'IID','MMN'}

        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = L*L';

end

% Inverse solution
%==========================================================================

% Greedy search over MSPs
%--------------------------------------------------------------------------
if strcmp(type,'GS')
    
    % extract patterns
    %----------------------------------------------------------------------
    Np    = length(Qp);
    Q     = sparse(Ns,Np);
    for i = 1:Np
        Q(:,i) = Qp{i}.q;
    end

    % Multivariate Bayes
    %----------------------------------------------------------------------
    MVB   = spm_mvb(UY,L,[],Q,Qe,16,1/2);
    M     = MVB.M;
    Cq    = MVB.qC;
    F     = max(MVB.F);

else
    
    % or ReML - ARD
    %----------------------------------------------------------------------
    Q           = {Qe{:} LQpL{:}};
    [Cy,h,Ph,F] = spm_sp_reml(YY,[],Q,Nr*Nt);

    % Covariances: sensor space - Ce and source space - L*Cp
    %----------------------------------------------------------------------
    Ne    = length(Qe);
    Np    = length(Qp);
    Ce    = sparse(Nm,Nm);
    LCp   = sparse(Nm,Ns);
    he    = h([1:Ne]);
    hp    = h([1:Np] + Ne);
    for i = 1:Ne
        try
            Ce = Ce + he(i)*Qe{i}.q*Qe{i}.q';
        catch
            Ce = Ce + he(i)*Qe{i};
        end
    end
    for i = 1:Np
        try
            LCp = LCp + hp(i)*LQpL{i}.q*Qp{i}.q';
        catch
            LCp = LCp + hp(i)*L*Qp{i};
        end
    end

    % MAP estimates of instantaneous sources
    %======================================================================
    iC    = inv(Cy);
    M     = LCp'*iC;

    % conditional covariance (leading diagonal)
    fprintf('Computing conditional covariances:')
    % Cq    = Cp - Cp*L'*iC*L*Cp;
    %----------------------------------------------------------------------
    Cp    = sparse(Ns,1);
    for i = 1:Np
        try
            Cp = Cp + hp(i)*sum(Qp{i}.q.^2,2);
        catch
            Cp = Cp + hp(i)*diag(Qp{i});
        end
    end
    Cq    = Cp - sum(LCp.*M')';

end

% select the most 'energetic' dipoles and store their indices in j
%--------------------------------------------------------------------------
[i j] = sort(-sum(M.^2,2));
j     = j(1:Na);
M     = M(j,:);
L     = L(:,j);
Is    = Is(j);
Cq    = Cq(j);

fprintf(' - done\n')

% re-scale and evaluate conditional expectation (of the sum over trials)
%--------------------------------------------------------------------------
SSR   = 0;
SST   = 0;
for i = 1:Nt
    
    Y{i}  = U'*Y{i};
    J{i}  = M*Y{i};

    % sum of squares
    %----------------------------------------------------------------------
    SSR   = SSR + sum(var((Y{i} - L*J{i}),0,2));
    SST   = SST + sum(var(Y{i},0,2));

end

% Assess accuracy; signal to noise (over sources)
%==========================================================================
R2    = 100*(SST - SSR)/SST;
fprintf('Variance explained %.2f (percent)\n',R2)

% Save results
%==========================================================================
model.inverse.type   = type;                 % inverse model
model.inverse.smooth = s;                    % smoothness (0 - 1)
model.inverse.xyz    = xyz;                  % VOI (XYZ)
model.inverse.rad    = rad;                  % VOI (rad)

model.inverse.M      = M;                    % MAP projector (reduced)
model.inverse.J      = J;                    % Conditional expectation
model.inverse.Y      = Y;                    % ERP data (reduced)
model.inverse.L      = L;                    % Lead-field  (reduced)
model.inverse.R      = speye(Nc,Nc);         % Re-referencing matrix
model.inverse.qC     = Cq;                   % spatial  covariance
model.inverse.qV     = Vq;                   % temporal correlations
model.inverse.T      = T;                    % temporal subspace
model.inverse.U      = U;                    % spatial  subspace
model.inverse.Is     = Is;                   % Indices of active dipoles
model.inverse.It     = It;                   % Indices of time bins
model.inverse.Ic     = Ic;                   % Indices of good channels
model.inverse.Nd     = Nd;                   % number of dipoles
model.inverse.pst    = pst;                  % pers-stimulus time
model.inverse.dct    = dct;                  % frequency range
model.inverse.F      = F;                    % log-evidence
model.inverse.R2     = R2;                   % variance accounted for (%)
model.inverse.woi    = woi;		     % timewindow inverted




% save in struct
%--------------------------------------------------------------------------
D.inv{D.val}         = model;
D.inv{D.val}.method  = 'Imaging';

% and delete old contrasts
%--------------------------------------------------------------------------
try
    D.inv{D.val} = rmfield(D.inv{D.val},'contrast');
end


% display
%==========================================================================
spm_eeg_invert_display(D);

return





% NOTES
%==========================================================================

% informed spatial basis functions using a generlised eigemode solution
%--------------------------------------------------------------------------
Cpos  = Y*(T'*diag(pst > 0)*T)*Y' + speye(Nc,Nc)*exp(-8);
Cpre  = Y*(T'*diag(pst < 0)*T)*Y' + speye(Nc,Nc)*exp(-8);
[U S] = spm_svd(inv(Cpos)*(Cpre),0);
U     = U(:,end - Nm + 1:end);



% display a selected basis function
%--------------------------------------------------------------------------
ql    = QG(:,Ip(64));
il    = find(ql > 1/64);
subplot(2,1,1)
spm_mip(ql(il),vert(il,:)',6);
axis image
