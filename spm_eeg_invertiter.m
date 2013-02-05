function [Dtest,modelF,allF]=spm_eeg_invertiter(Dtest,Np,Npatchiter,Nm,Nt)

%  Function to perform several MSP type inversions with different
%  pseudo-randomly selected priors- in this case single cortical patches
% Np :number of patches to be used per iteration
% Npatchiter: number of iterations
% Nm: number of spatial modes 
% Nt: number of temporal modes
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
% 
% Gareth Barnes
% $Id: spm_eeg_invertiter.m 5244 2013-02-05 17:05:47Z gareth $


disp('running iterative inversion');
if nargin<4,
    disp('estimating spatial modes from lead fields');
    Nm=[];
end;

if nargin<5,
    disp('estimating temporal modes from data');
    Nt=[];
end;


if ~isfield(Dtest{1},'val'),
    val=1;
end;
    

Nvert=size(Dtest{1}.inv{val}.mesh.tess_mni.vert,1);
twindow=Dtest{1}.inv{val}.inverse.woi./1000; %% in sec
fband=[Dtest{1}.inv{val}.inverse.hpf Dtest{1}.inv{1}.inverse.lpf];
allF=zeros(Npatchiter,1);
disp('Reseting random number seed !');
rand('state',0); 
for patchiter=1:Npatchiter, %% change patches
    
    randind=randperm(Nvert);
    Ip=randind(1:Np);
    
    
    
    
    Dtest{1}	= spm_eeg_invert_noscale(Dtest{1},Ip,fband,Nm,Nt,[],twindow);			% Source reconstruction
    Dtest{1}.inv{val}.inverse.Ip=Ip;
    if isempty(Nm),
        allNt=max([size(Dtest{1}.inv{val}.inverse.T,2),1]);
        allNm=size(Dtest{1}.inv{val}.inverse.M,2);
        Nm=allNm;
        Nt=allNt;
    end;
    
    modelF(patchiter).inverse=Dtest{1}.inv{val}.inverse;
    allF(patchiter)=[Dtest{1}.inv{val}.inverse.F];
    allJ{patchiter}=Dtest{1}.inv{val}.inverse.J;
    manyinverse{patchiter}=modelF(patchiter).inverse;
    disp('iteration:');
    patchiter
    allF
    
    
end; % for patchiter

[bestF,bestind]=max(allF);
disp('model evidences relative to maximum:')

sort(allF-bestF)

Dtest{1}.inv{val}.inverse=modelF(bestind).inverse; %% return best model for now
if Dtest{1}.inv{val}.inverse.BMA,
    disp('Running BMA to get current estimate');
    Jbma=spm_bma_msp(manyinverse,allF); %% onlt the mean is calculated using BMA (but could extend to covariance)
    Dtest{1}.inv{val}.inverse.T=1; %% Jbma is the sum of all modes
    Dtest{1}.inv{val}.inverse.J={Jbma};
else
    disp('Using best patch set to current estimate');
end; % if BMA

Dtest{1}.inv{val}.inverse.allF=allF;
spm_eeg_invert_display(Dtest{1});
%Dtest{1}.inv{val}.inverse.J=Jbma;



  




