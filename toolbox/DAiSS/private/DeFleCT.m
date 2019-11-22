function w=DeFleCT(passband,SVDpassband,force_passband_flag,stopband,SVDstopband,...
    LFM,C,SNR2,Csvdtrunc)
% function w=DeFleCT(passband,SVDpassband,force_passband_flag,stopband,SVDstopband,...
%     LFM,C,SNR2,Csvdtrunc)
% function w=DeFleCT(passband,SVDpassband,force_passband_flag,stopband,SVDstopband,...
%     LFM,[],SNR2,Whitener)
%
% Makes a DeFleCT spatial filter with given passband and stopband. 
% passband:     indices to sources for which the targeted output is 1
% SVDpassband:  how many components represent the passband (optional)
% force_passband_flag:  forces the output for all passband components to 1
% stopband:             indices to sources for which the output is 0
% SVDstopband:  how many components represent the stopband (optional)
% LFM:  forward model
% C:    noise covariance matrix (or measurement covariance matrix)
% SNR2: assumed signal-to-noise ratio (for setting regularization)
% Csvdtrunc: number of truncated singular values for the inversion of C.
% Whitener:  the whitener that is applied to the leadfields (optional; if C
%            is given, the whitener is built from C)
%
% L, C, SNR2, Csvdtrunc/Whitener need to be given only, when any of these
% parameters changes
%
% This function is implemented according to:
% Hauk and Stenroos: A Framework for the Design of Flexible Cross-Talk
% Functions for Spatial Filtering of EEG/MEG Data: DeFleCT. HBM 2013.
%
% Copyright (c) 2012--2013 Matti Stenroos (matti.stenroos@aalto.fi)
% -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
%                !! There is no warranty of any kind !!
% -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
% 17 Oct 2013


persistent Sinv L Nch Ns W
if nargin>5
    [Nch,Ns]=size(LFM);
    if length(Csvdtrunc)==1
        W=MakeWhitener(C,Csvdtrunc);
    else
        W=Csvdtrunc;
    end
    L=W*LFM;%persistent
    lambda2=trace(L*L')/(Nch*SNR2);
    S=L*L'+lambda2*eye(Nch);
    Sinv=inv(S);%persistent
end


%how many indices in the stopband and passband
if nargin>3, Nstop=length(stopband);end
Npass=length(passband);

if nargin<3 || isempty(force_passband_flag), force_passband_flag=0;end

%Build projector P, wP=i strictly.
%if force_passband_flag, then set wP=1 for the passband vector, else leave
%passband out.
Ppass=L(:,passband);
Pstop=L(:,stopband);

%take svd components for the stopband
if nargin>4 && ~isempty(SVDstopband) && SVDstopband<Nstop
    [U,s,Vt]=svd(Pstop,0);
    Pstop=U(:,1:SVDstopband)*s(1:SVDstopband,1:SVDstopband);
    Nstop=SVDstopband;
end
%take svd components for the passband
if nargin>1 && ~isempty(SVDpassband) && SVDpassband<Npass
    [U,s,Vt]=svd(Ppass,0);
    Ppass=U(:,1:SVDpassband)*s(1:SVDpassband,1:SVDpassband);
    Npass=SVDpassband;
end
%make strict projector, with or without passband
if nargin>2 && force_passband_flag
    P=[Ppass Pstop];
    i=zeros(1,Npass+Nstop);
    i(1:Npass)=1;
else
    P=Pstop;
    i=zeros(1,Nstop);
end
%make passband-vector for "traditional" minimisation
t=zeros(1,Ns);
if ~force_passband_flag
    t(passband)=1;
end
% and make the final minimiser
% Pinv=inv(P'*Sinv*P);
% term2=(i-t*L'*Sinv*P)*Pinv*P';
term2=(i-t*L'*Sinv*P)/(P'*Sinv*P)*P';
w=(t*L'+term2)*Sinv*W;

function [W,deW]=MakeWhitener(C,trunc)
% function [W,deW]=MakeWhitener(C,trunc)
% Computes whitener C^(-1/2) and de-whitener sqrt(C)
% trunc: how many smallest eigenvalues are set to zero before inversion
% (all other regularization has to be applied to C before calling this
% routine)
%
% Copyright (c) 2012--2013 Matti Stenroos (matti.stenroos@aalto.fi)
% -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
%                !! There is no warranty of any kind !!
% -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
% version 130610
if min(size(C))==1,
    C=diag(diag(C));
end

[U,D]=eig(C);
[ds,sortind]=sort(diag(D));
if trunc
    ds(1:trunc)=0;
end
Us=U(:,sortind);
dsinv=zeros(size(ds));
dsinv((trunc+1):end)=1./ds((trunc+1):end);
W=Us*diag(sqrt(dsinv))*Us';

if nargout==2
    deW=U*diag(sqrt(diag(D)))*Us';
end
