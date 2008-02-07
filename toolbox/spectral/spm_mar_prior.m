function [prior] = spm_mar_prior (d,p,type)
% Specify ARD-type prior for Bayesian MAR model
% function [prior] = spm_mar_prior (d,p,type)
%
% d         Number of time series
% p         Order of MAR model
% type      'global', 'lag','interaction','lag-inter',
%           'silly','ran2','triu' (see code below)
%
% prior     data structure to be passed to spm_mar.m
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_mar_prior.m 1143 2008-02-07 19:33:33Z spm $

k=p*d*d;
switch type,
    case 'global',
        prior.type='global';
        prior.groups=1;
        prior.group(1).index=ones(k,1);
    case 'ran2',
        % Randomly group together weights into 2 groups
        k1=floor(k/2);
        prior.type='ran2';
        prior.groups=2;
        Ij=randperm(k);
        prior.group(1).index=zeros(k,1);
        prior.group(1).index(Ij(1:k1))=ones(1,k/2);
        prior.group(2).index=zeros(k,1);
        prior.group(2).index(Ij(k1+1:k))=ones(1,k/2);
    case 'silly',
        % Randomly group together 2 weights
        prior.type='silly';
        prior.groups=k/2;
        Ij=randperm(k);
        for j=1:k/2,
            prior.group(j).index=zeros(k,1);
            prior.group(j).index(Ij(j))=1;
            prior.group(j).index(Ij(j+k/2))=1;
        end
    case 'lag',
        % p groups; one for each lag
        prior.type='lag';
        prior.groups=p;
        for j=1:p,
            Iwj=zeros(p*d,d);
            start=(j-1)*d+1;
            stop=start+d-1;
            Iwj(start:stop,:)=ones(d,d);
            prior.group(j).index=Iwj(:);
        end
    case 'interaction',
        % (1) One group for within-series coefficients
        % (2) Another group for between-series coefficients
        prior.type='interaction';
        prior.groups=2;
        Iw=[];
        for i=1:p,
            Iw=[Iw;diag(ones(1,d))];
        end
        prior.group(1).index=Iw(:);
        prior.group(2).index=ones(k,1)-Iw(:);
    case 'triu',
        % (1) One group for upper-triangular coefficients
        %     - including diagonal
        % (2) Another group for lower
        prior.type='triu';
        prior.groups=2;
        Iw=[];
        for i=1:p,
            Iw=[Iw;triu(ones(d,d))];
        end
        prior.group(1).index=Iw(:);
        prior.group(2).index=ones(k,1)-Iw(:);
    case 'tril',
        % (1) One group for lower-triangular coefficients
        %     - including diagonal
        % (2) Another group for upper
        prior.type='tril';
        prior.groups=2;
        Iw=[];
        for i=1:p,
            Iw=[Iw;tril(ones(d,d))];
        end
        prior.group(1).index=Iw(:);
        prior.group(2).index=ones(k,1)-Iw(:);
    case 'lag-inter',
        % Lag-interaction prior
        % 2p groups.
        % At each lag have 
        % (a) a group for within series coeffs and 
        % (b) a group for interaction coeffs
        prior.type='lag-inter';
        prior.groups=2*p;
        % First p groups are the within terms
        % Next p groups are the interaction terms
        for j=1:p,
            Iwj_within=zeros(p*d,d);
            Iwj_inter=zeros(p*d,d);
            start=(j-1)*d+1;
            stop=start+d-1;
            Iwj_within(start:stop,:)=diag(ones(1,d));
            Iwj_inter(start:stop,:)=ones(d,d)-diag(ones(1,d));
            prior.group(j).index=Iwj_within(:);
            prior.group(j+p).index=Iwj_inter(:);
        end
        
    otherwise,
        disp('Unknown prior type in marprior');
        prior=[];
        return
end
