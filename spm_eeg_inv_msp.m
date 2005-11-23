function [varargout] = spm_eeg_inv_msp(varargin)

%=======================================================================
% Multivariate Source Prelocalization
% see Mattout et al., NeuroImage 2005 Jun 26 2:356-73.
%
% FORMAT [APM,Cvel,Nfact] = spm_eeg_inv_msp(Gnorm,VectP,ValP)
% Input:
% Gnorm		- normalized gain matrix with normalized
% VectP     - matrix containing the eigenvectors of Gnorm
% ValP      - diagonal matrix containing the eigenvalues of Gnorm
% (see spm_eeg_inv_PCAgain.m)
% Output:
% APM       - vector of probability of activation (size = #dipoles)
% Cvel      - Velicer criterion
% Nfact     - #eigenvectors selected by Velicer filtering
% Ce        - estimated noise covariance matrix by projecting the data
%             onto Velicer's subspace
%
% FORMAT [fname_out] = spm_eeg_inv_msp(fname_in)
% Input:
% fname_in  - .mat file containing the output of spm_eeg_inv_PCAgain.m
% Output:
% fname_out - .mat file containing the output variables
%
% FORMAT D = spm_eeg_inv_msp(S)
% Input:
% S     - EEG/MEG data structure
% Output:
% D     - same data structure including the new output parameters/files
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_msp.m 308 2005-11-23 19:21:56Z jeremie $

spm_defaults

epsilon = 1e10*eps;

if nargin == 1
    if ischar(varargin{1})
        variabl = load(varargin{1});
        names   = fieldnames(variabl);
        Gnorm   = getfield(variabl,'Gnorm');
        VectP   = getfield(variabl,'VectP');
        ValP    = getfield(variabl,'ValP');
        clear variabl
    else
        try
            D = varargin{1};
        catch
            D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
        	D = spm_eeg_ldata(D);
        end
        val = length(D.inv);
        variabl = load(D.inv{val}.forward.pcagain);
        names   = fieldnames(variabl);
        Gnorm   = getfield(variabl,'Gnorm');
        VectP   = getfield(variabl,'VectP');
        ValP    = getfield(variabl,'ValP');
        clear variabl
    end
elseif nargin == 3
    Gnorm = varargin{1};
    VectP = varargin{2};
    ValP  = varargin{3};
else    
    disp('Wrong input format');
    return
end
    
woi       = D.inv{val}.inverse.woi;
contrast  = D.inv{val}.inverse.contrast;
Nsens     = D.Nchannels;

woi(1) = round(woi(1)*(D.Radc/1000)) + D.events.start + 1;
woi(2) = round(woi(2)*(D.Radc/1000)) + D.events.start + 1;
Y      = sparse(zeros(size(D.data,1),(woi(2)-woi(1)+1)));;
for i = 1:D.Nevents
    k = find(D.events.code == D.events.types(i));
    Y = Y + contrast(i)*squeeze(D.data(:,woi(1):woi(2),k));
end

% Data normalization
NCM   = norcol(Y);
O     = ones(Nsens,1);
Lm    = O*NCM;
Ynorm = Y./Lm;
Stot  = sum(ValP.^2);
clear ValP
 
% Velicer criterion based filtering
[Nfact,Cvel] = spm_eeg_inv_velicer(VectP, Ynorm);

Bs      = VectP(:,1:Nsens-Nfact);
Gamma_s = Bs'*Ynorm;
Ys      = Bs*Gamma_s;
Op      = pinv(Gamma_s'*Gamma_s);
Op      = Ys*Op*Ys';
Gs      = Op*Gnorm;
APM     = sum(Gs.*Gs);
APM     = epsilon + (1 - 2*epsilon)*(APM-min(APM))/(max(APM)-min(APM));
APM     = APM.^2;

% Estimating the noise covariance matrix
SubS    = VectP(:,Nfact+1:end);
Gamma_n = SubS'*Y;
Yn      = SubS*Gamma_n;
Ce      = Yn*Yn'/(size(Yn,2));

if nargin == 1
    woi = D.inv{val}.inverse.woi;
    if ischar(varargin{1})
        [pth,nam,ext] = spm_fileparts(varargin{1});
        Idelim    = find(nam == '_');
        fname_out = nam(1:Idelim(end-1)-1);
        fname_out = [fname_out '_msp_' num2str(woi(1)) '-' num2str(woi(2)) 'ms.mat'];
        save(fullfile(pth,fname_out),'APM','Cvel','Nfact','Ce');
        varargout{1} = fname_out;
    else
        [pth,nam,ext] = spm_fileparts(D.inv{val}.forward.pcagain);
        Idelim    = find(nam == '_');
        fname_out = nam(1:Idelim(end-1)-1);
        fname_out = [fname_out '_msp_' num2str(woi(1)) '-' num2str(woi(2)) 'ms.mat'];
        D.inv{val}.inverse.priors.level2{3}.filename = fullfile(pth,fname_out);
        save(D.inv{val}.inverse.priors.level2{3}.filename,'APM','Cvel','Nfact','Ce');
        save(fullfile(D.path, D.fname), 'D');
        varargout{1} = D;
    end
else
    if nargout == 4
        varargout{1} = APM;
        varargout{2} = Cvel;
        varargout{3} = Nfact;
        varargout(4) = Ce;
    else
        disp('Wrong output format');
    end
end


%=======================================================================
function [Nfact,Cvel] = spm_eeg_inv_velicer(B,Y)
% The velicer criterion enables to estimate the data subspace whose
% covariance is the closest to a diagonal matrix (uncorrelated noise).
% This subspace is iteratively built with the eigenvectors of the data
% covarance matrix which correspond to the smallest eigenvalues.
%
% see [Velicer W.F., 1976, Determining the number of components from matrix
% of partial correlations. Psychometrika, 41:321-327]

dimt  = size(Y,1);
tim   = size(Y,2);

Cvel = zeros(1,dimt-1);
for i = 1:(dimt-1)
    Gam = B(:,dimt-i:dimt)'*Y;
    Cov = Gam'*Gam;
    Cov = Cov.^2;
    d = diag(Cov);
    Y1 = d(:,ones(tim,1));
    Y2 = Y1';
    Y3 = sqrt(Y1 .* Y2);
    Cov = Cov ./ Y3;
    for k = 1:tim
        Cov(k,k) = 0;
    end
    Cvel(i) = sum(sum(Cov))/(dimt*(dimt-1));
end

Nfact = find(Cvel == min(Cvel));

