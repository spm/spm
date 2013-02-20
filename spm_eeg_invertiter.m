function [Dtest,modelF,allF]=spm_eeg_invertiter(Dtest,Npatchiter,funcname)

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
% $Id: spm_eeg_invertiter.m 5268 2013-02-20 14:46:38Z gareth $

if nargin<2,
    Npatchiter=[];
end;


if nargin<3,
    error('No function call specified');
end;

if isempty(Npatchiter),
    Npatchiter=16;
end;

if numel(Dtest)>1,
    error('only works with single datasets at the moment');e
end;

val=Dtest{1}.val;
Nvert=size(Dtest{1}.inv{val}.mesh.tess_mni.vert,1);
Np=Dtest{1}.inv{val}.inverse.Np;

allF=zeros(Npatchiter,1);
disp('Reseting random number seed !');
rand('state',0);
for patchiter=1:Npatchiter, %% change patches
    
    randind=randperm(Nvert);
    Ip=randind(1:Np);
    
    
    switch funcname,
        case 'Classic',
            
            Dtest{1}.inv{val}.inverse.Ip=Ip;
            
            Dtest{1}	= spm_eeg_invert_classic(Dtest{1});
        case 'Current'
            warning('Patch centres are currently fixed for this algorithm (iteration will have no effect!)');
            
    
            Dtest{1}	= spm_eeg_invert(Dtest{1}); %
            Dtest{1}.inv{val}.inverse.Ip=Ip;
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
    [Jbma,qCbma]=spm_eeg_invert_bma(manyinverse,allF); %% onlt the mean is calculated using BMA (but could extend to covariance)
    %Dtest{1}.inv{val}.inverse.T=1; %% Jbma is the sum of all modes
    Dtest{1}.inv{val}.inverse.J={Jbma};
    Dtest{1}.inv{val}.inverse.qC=qCbma;
    
else
    disp('Using best patch set to current estimate');
end; % if BMA

keyboard
Dtest{1}.inv{val}.inverse.BMAF=allF;
spm_eeg_invert_display(Dtest{1});









