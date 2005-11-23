function D = spm_eeg_inv_inverse(S,varargin)

%=======================================================================
% Compute the requested variance components at both sensor and source level
%
% FORMAT D = spm_eeg_inv_inverse(S)
% Input:
% S		    - input data struct (optional)
% Output:
% D			- same data struct including the inverse solution files and variables
%
% FORMAT D = spm_eeg_inv_priors(S,Vsens,Vsour)
% Input:
% S		    - input data struct (optional)
% Vsens     - vector indicating the number of component to use at the
%             sensor level (default [1 0 0]).
%             Vsens(1) relates to the iid component
%              (=1 component is included, =0 it is not)
%             Vsens(2) relates to an anti-average estimate of the noise
%              (only if evoked response)
%             Vsens(3) relates to an empirical estimate/filtering of the
%             data (only if evoked response)
%             Additive component are specified by the user
% Vsour     - vector indicating the number of priors to use at the source
%             level (default [1 0 0]).
%             Vsour(1) relates to a smoothness constraint
%              (=1 it is taken into account, =0 it is not)
%             Vsour(2) relates to a minimum norm constraint
%             Vsour(3) relates to a Multivariate Source Prelocalisation
%             constraint (only if evoked response)
%             Additive priors can be specified by the user
% Output:
% D			- same data struct including the inverse solution files and variables
%
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id$

spm_defaults
s = 8; % smoothness parameter (millimeters)

try
    D = S; 
    clear S
catch
    D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
	D = spm_eeg_ldata(D);
end

if nargin == 3
    Vsens = varargin{1};
    Vsour = varargin{2};
elseif nargin == 1
    Vsens = [1 0 0];
    Vsour = [1 0 0];
else
    error(sprintf('Wrong input arguments\n'));
end  

val = length(D.inv);
Nv  = D.inv{val}.mesh.Ctx_Nv;
Ns  = D.Nchannels;

% LOAD/COMPUTE NOISE COVARIANCE COMPONENTS - SENSOR SPACE (level 1)
Qe    = {};

if Vsens(1)     % i.i.d component
    Qe{1} = speye(Ns);
    D.inv{val}.inverse.priors.level1{1}.status = 'yes';
else
    D.inv{val}.inverse.priors.level1{1}.status = 'no';
end

if Vsens(2)     % Anti-average estimate of noise covariance matrix
    D.inv{val}.inverse.priors.level1{2}.status = 'yes';
    if isempty(D.inv{val}.inverse.priors.level1{2}.filename)
        Qe{Vsens(1) + 1} = spm_eeg_inv_AntiAverage(D.fname,D.inv{val}.inverse.contrast,D.inv{val}.inverse.woi);
    else
        variabl = load(D.inv{val}.inverse.priors.level1{2}.filename);
        name = fieldnames(variabl); 
        Qe{Vsens(1) + 1} = getfield(variabl , name{1});
    end
else
    D.inv{val}.inverse.priors.level1{2}.status = 'no';
end

Ncomp = length(Vsens) - 3;
if Ncomp > 0
    np = sum(Vsens(1:3));
    for i = 1:Ncomp
        if Vsens(3+i)
            if ~isempty(D.inv{val}.inverse.priors.level1{3+i}.filename)
                variabl = load(D.inv{val}.inverse.priors.level1{3+i}.filename);
                name = fieldnames(variabl);
                if isempty(D.inv{val}.inverse.priors.level1{3+i}.label)
                    D.inv{val}.inverse.priors.level1{3+i}.label = name{1};    
                end
                Qe{np+i} = getfield(variabl , name{1});
                D.inv{val}.inverse.priors.level1{3+i}.status = 'yes';
                if any(size(Qe{np+i}) - [Ns Ns])
                    error(sprintf('Wrong matrix size\n'));
                end
            else
                error(sprintf('No file to download\n'));
            end
        else
            D.inv{val}.inverse.priors.level1{3+i} = 'no';
        end
    end
end


% LOAD/COMPUTE PRIOR VARIANCE COMPONENTS - SOURCE SPACE (level 2)
Qp     = {};

if Vsour(1) % Compute Smoothness Prior
    if ~isempty(D.inv{val}.mesh.CtxGeoDist)
        variabl = load(D.inv{val}.mesh.CtxGeoDist);
        name    = fieldnames(variabl);
        Mdist   = getfield(variabl , name{1});
    else
        try
            load(D.inv{val}.mesh.tess_ctx);
            Mdist = spm_eeg_inv_meshdist(vert,face);
            clear vert face
        catch
            disp('Missing Mesh file');
            return
        end
    end
    if any(size(Mdist) - [Nv Nv])
        disp('Wrong mesh size');
        return
    end
    Qp{1} = speye(size(Mdist));
    Mdist = Mdist.^2/(2*s^2);
    [I,J] = find(Mdist);
    K = I - J;
    if ~all(K)
        IndZ = find(K == 0);
        I(IndZ) = [];
        J(IndZ) = [];
    end
    clear K
    for k = 1:length(I)
            Qp{1}(I(k),J(k)) = exp(-Mdist(I(k),J(k)));
    end        
    clear Mdist
    D.inv{val}.inverse.priors.level2{1}.status = 'yes';
else
    D.inv{val}.inverse.priors.level2{1}.status = 'no';
end

if Vsour(2) % Minimum Norm prior
    Qp{Vsour(1) + 1} = speye(Nv);
    D.inv{val}.inverse.priors.level2{2}.status = 'yes';
else
    D.inv{val}.inverse.priors.level2{1}.status = 'no';
end

if Vsour(3) | Vsens(3) | (Nv - D.inv{val}.inverse.dim)  % Mutlivariate Source Prelocalisation (MSP) prior
    if ~isempty(D.inv{val}.inverse.priors.level2{3}.filename)
        variabl = load(D.inv{val}.inverse.priors.level2{3}.filename);
        name = fieldnames(variabl);
        if Vsens(3)
            Ce = getfield(variabl , 'Ce');
            Qe{Vsens(1) + Vsens(2) + 1} = Ce;
            D.inv{val}.inverse.priors.level1{3}.status = 'yes';
        else
            D.inv{val}.inverse.priors.level1{3}.status = 'no';
        end
        if Vsour(3)
            APM = getfield(variabl , 'APM');
            Qp{Vsour(1) + Vsour(2) + 1} = diag(APM);
            D.inv{val}.inverse.priors.level2{3}.status = 'yes';
            clear APM
        else
            D.inv{val}.inverse.priors.level2{3}.status = 'no';
        end
    else    
        if ~isempty(D.inv{val}.forward.pcagain)
            D = spm_eeg_inv_msp(D);
        elseif ~isempty(D.inv{val}.forward.gainmat)
            variabl = load(D.inv{val}.forward.gainmat);
            name    = fieldnames(variabl);
            G       = getfield(variabl , name{1});
            [Gnorm,VectP,ValP] = spm_eeg_inv_PCAgain(G);
            [pth,nam,ext] = spm_fileparts(variabl);
            D.inv{val}.forward.pcagain = fullfile(pth,[nam '_pca.mat']);
            save(D.inv{val}.forward.pcagain,'Gnorm','VectP','ValP');
            D = spm_eeg_inv_msp(D);
            clear variabl Gnorm VectP ValP G
        else
            error(sprintf('No gain matrix available\n'));
        end
    
        variabl = load(D.inv{val}.inverse.priors.level2{3}.filename);
        name    = fieldnames(variabl);
    
        if Vsens(3)
            Ce  = getfield(variabl , 'Ce');
            Qe{Vsens(1) + Vsens(2) + 1} = Ce;
            D.inv{val}.inverse.priors.level1{3}.status = 'yes';
        else
            D.inv{val}.inverse.priors.level1{3}.status = 'no';
        end
        
        if Vsour(3)
            APM = getfield(variabl , 'APM');
            Qp{Vsour(1) + Vsour(2) + 1} = diag(APM);
            D.inv{val}.inverse.priors.level2{3}.status = 'yes';
        else
            D.inv{val}.inverse.priors.level2{3}.status = 'no';
        end
    
        clear variabl           
    end
end

Nprior = length(Vsour) - 3;
if Nprior > 0
    np = sum(Vsour(1:3));    
    for i = 1:Nprior
        if Vsour(3+i)
            if ~isempty(D.inv{val}.inverse.priors.level2{3+i}.filename)
                variabl = load(D.inv{val}.inverse.priors.level2{3+i}.filename);
                name = fieldnames(variabl);
                D.inv{val}.inverse.priors.level2{3+i}.label = name{1};
                Qp{np+i} = getfield(variabl , name{1});
                D.inv{val}.inverse.priors.level2{3+i}.status = 'yes';
                if any(size(Qp{np+i}) - [Nv Nv])
                    disp('Wrong matrix size');
                    break
                end
                clear variabl
            else
                error(sprintf('No file to download\n'));
            end
        else
            D.inv{val}.inverse.priors.level2{3+i}.status = 'no';
        end
    end
end


% INVERSE COMPUTATION
if D.inv{val}.inverse.activity == 'evoked'
    D = spm_eeg_inv_evoked(D,Qe,Qp);
elseif (D.inv{val}.inverse.activity == 'induced') | (D.inv{val}.inverse.activity == 'evoked & induced')
    D = spm_eeg_inv_induced(D,Qe,Qp)
else
    error(sprintf('Missing analysis type\n'));
end


%=======================================================================
function Mdist = spm_eeg_inv_meshdist(vert,face)
% Efficient computation of the 2nd order distance matrix of a triangulated
% irregular mesh, based on the cortical neighbourhoud (geodesic distance)
%
% Inspired by function mesh_laplacian.m by Darren Weber
% from the bioelectromagnetism matlab toolbox
% see http://eeg.sourceforge.net/

Nv = length(vert);
Nf = length(face);

edge = zeros(N);
for i = 1:size(face,1);
    Diff  = [vert(face(i,[1 2 3]),:) - vert(face(i,[2 3 1]),:)];
    EuclD = sqrt( sum(Diff.^2, 2) );
    
    edge(face(i,1),face(i,2)) = EuclD(1);
    edge(face(i,2),face(i,3)) = EuclD(2);
    edge(face(i,3),face(i,1)) = EuclD(3);
    
    edge(face(i,2),face(i,1)) = EuclD(1);
    edge(face(i,3),face(i,2)) = EuclD(2);
    edge(face(i,1),face(i,3)) = EuclD(3);
end
clear face vert

Mdist = edge;
for i = 1:Nv
    a          = find(edge(i,:));
    [b,c]      = find(edge(a,:));
    Mdist(i,c) = Mdist(i,a(b)) + diag(Mdist(a(b),c))';
    Mdist(c,i) = Mdist(i,c)';
end
clear edge

return
%=======================================================================


%=======================================================================
function Caa = spm_eeg_inv_AntiAverage(fname,contrast,woi)
% Compute an empirical estimate of the noise covariance matrix
% by anti-averaging the epoched data, thus removing the evoked component
% and keeping the random effect.

[pth,nam,ext] = spm_fileparts(fname);
Idelim        = find(nam == '_');
firstnam      = nam(1:Idelim(1)-1);
Im            = find(firstnam == 'm');
reducenam     = firstnam;
reducenam(Im) = [];
newnam        = [reducenam '_' nam(Idelim(1)+1:end) '.mat'];
fnewname      = fullfile(pth,fnewname);

try
    D         = spm_eeg_ldata(fnewname);
catch
    D         = spm_select(1, '.mat', 'Select non-average data file');
	D         = spm_eeg_ldata(D);
end

woi(1) = round(woi(1)*(D.Radc/1000)) + D.events.start + 1;
woi(2) = round(woi(2)*(D.Radc/1000)) + D.events.start + 1;

Inz    = find(contrast);
Ye     = [];
N      = 0;
T      = woi(2) - woi(1) + 1;

for i = 1:length(Inz)
    Trials = find(D.events.code == D.events.type(Inz(i)));
    n(i)   = length(Trials);
    if round(n(i)/2) ~= n(i)/2
        n(i) = n(i) - 1;
    end
    randvec = randperm(n(i));
    Yi     = [];
    for k = 1:n(i)
        s = 1;
        if round(k/2) == k/2
            s = -1;
        end
        Yi = [Yi s*squeeze(D.data(:,woi(1):woi(2),Trials(randvec(k))))];
    end
    Ye = Ye + Yi*kron(ones(n,1),speye(T));
    N  = N + n(i)
    clear Yi
end
clear D

Ye  = Ye/N;
Caa = Ye*Ye'/T;

return
%=======================================================================

