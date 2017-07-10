function [LCpL,Q,sumLCpL,QE,Cy,M,Cp,Cq,Lq] = spm_eeg_assemble_priors(L,Qp,Qe,ploton,h)
% Predict sensor level impact of sources in Qp given sensor noise Qe
% FORMAT [LCpL,Q,sumLCpL,QE,Cy,M,Cp,Cq,Lq] = spm_eeg_assemble_priors(L,Qp,Qe,ploton,h)
% L       - lead fields
% Qp      - priors on source level dipole moment nAm/(mm2) per sample
% Qe      - sensor level variance in fT^2 per sample
% h       - optional hyperparameters that scale the variance components in
%           Qe and Qp (assume sensor followed by source level parameters)

% LCpL    - sensor level covariance components corresponding to the source
%           priors Qp
% Q       - sparse array over sources holding dipole moment density nAm/(mm2)
% sumLCpL - predicted sensor level variance due to sources (in fT^2)
% QE      - predicted sensor level noise variance (in fT^2)
% Cy      - total sensor noise covariance predicted : Cy = QE+sumLCpL
% M       - MAP estimator : M = Cp*L'*inv(Qe + L*Cp*L'))
% Cp      - total source covariance matrix
% Cq      - conditional source covariance - need to implement
% Lq      - cell array of the L*q (impact of each source component at sensor level)
%__________________________________________________________________________
% Copyright (C) 2010-2017 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_eeg_assemble_priors.m 7132 2017-07-10 16:22:58Z guillaume $

if nargin < 4
    ploton = ~spm('CmdLine');
end

if nargin < 5
    h=[];
end

Ne = length(Qe);
if iscell(Qp)
    Np = length(Qp); % number of source components
else
    Np = size(Qp,2);
end

if isempty(h)
    h = ones(Ne+Np,1);
end

if length(h)~=Ne+Np
    error('There must be as many hyperparameters as priors.');
end

he = h(1:Ne); % sensor level
hp = h(Ne + (1:Np)); % source level

Ns = size(L,2); % number sources

Q = sparse(Ns,Np);
Lq = {};

fprintf('Assembling %d prior components\n',Np);
for i = 1:Np
    if iscell(Qp)
        if isfield(Qp{i},'q')
            %Q(:,i) = Qp{i}.q; % patch amplitudes can be positive or negative
            q = Qp{i}.q;
            v = 1;
            if isfield(Qp{i},'v')
                v = sqrt(Qp{i}.v);
            end
            
            Q(:,i) = q*v;
            
            % Lq is the sensor level projection of the prior Q{i}.q
            Lq{i}.q = L*Q(:,i); %  supply an eigen mode in q
        else % no .q field, check it is 2D
            if size(Qp{i},1)~=size(Qp{i},2)
                disp('making Qp 2D');
                Qp{i} = diag(sparse(Qp{i}));
            end
        end 
    end
end


%-Assemble empirical priors
%==========================================================================

LCpL  = {};

sumLCpL = zeros(size(L,1),size(L,1));

Cp  = sparse(0);
LCp = sparse(0);

for i = 1:Np
    
    if isfield(Qp{i},'q') % sparse diagonbal representation
        
        Q(:,i) = Q(:,i)*sqrt(hp(i)); % update moment by sqrt of variance scaling
        
        LQp = L*Q(:,i);
        LCpL{i} = LQp*LQp';
        Cp  = Cp + Q(:,i)*Q(:,i)'; %% Q is already scaled by root hp
        %LCp = LCp+LQp;
        LCp = LCp+L*Cp;
    else % full matrix version
        Qtmp = Qp{i}*hp(i);
        LCpL{i} = L*Qtmp*L';
        Cp  = Cp + Qtmp;
       % Q(:,i)=diag(Qtmp);
    end
    
    sumLCpL = sumLCpL + LCpL{i};
    
end

QE = zeros(size(Qe{1}));
for i=1:length(Qe)
    QE = QE + Qe{i}*he(i);
end

Cy = sumLCpL + QE;

% This is equivalent to M = Cp*UL'*inv(Qe + UL*Cp*UL'))
% with Cp the posterior source covariance (with optimal h values)
%Qsum=sum(Q(:,i).^2); % diagonal of posterior source cov matrix
%Cp=diag(Qsum);
M = Cp*L'/Cy;

% conditional variance (leading diagonal)
% Cq    = Cp - Cp*L'*iC*L*Cp;
%----------------------------------------------------------------------
%Cq    = Cp - sum(LCp.*M')'; in original

Cq = diag(Cp) - sum(LCp.*M')';

disp('end assemble');

if ploton
    subplot(3,1,1);
    imagesc(sumLCpL);colorbar;
    title('source prior');
    subplot(3,1,2);
    imagesc(QE);colorbar;
    title('sensor noise prior');
    subplot(3,1,3);
    imagesc(Cy);colorbar;
    title('Sensor covariance per sample (fT^2)');
end
