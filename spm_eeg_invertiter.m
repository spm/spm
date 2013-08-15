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
% $Id: spm_eeg_invertiter.m 5615 2013-08-15 14:37:24Z spm $

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
    error('only works with single datasets at the moment');
end;

val=Dtest{1}.val;
Nvert=size(Dtest{1}.inv{val}.mesh.tess_mni.vert,1);
Np=Dtest{1}.inv{val}.inverse.Np;

allF=zeros(Npatchiter,1);
disp('Reseting random number seed !');




fprintf('Checking leadfields');

try
    par = matlabpool('size') > 0;
catch
    par=0;
end;
if par,
    [L Dtest{1}] = parfor_spm_eeg_lgainmat(Dtest{1});
    
else
    [L Dtest{1}] = spm_eeg_lgainmat(Dtest{1});  % Generate/load lead field- this stops it being done at each iteration
end;


rand('state',0);


%% make random patch centers, except on the first iteration when we keep to a fixed stepsize
for f=1:Npatchiter,
    tmp=randperm(Nvert);
    allIp(f).Ip=tmp(1:Np);
end;
allIp(1).Ip=[]; %% fix the step size of the first patch set


if par, %% RUN IN PARALLEL
    parfor patchiter=1:Npatchiter,
        Din=Dtest{1};
        Dout=Din;
        Ip=allIp(patchiter).Ip;
        
        switch funcname,
            case 'Classic',
                Din.inv{val}.inverse.Ip=Ip;
                Dout    = spm_eeg_invert_classic(Din);
            case 'Current',
                warning('Patch centres are currently fixed for this algorithm (iteration will have no effect!)');
                Dout    = spm_eeg_invert(Din); %
                Dout.inv{val}.inverse.Ip=Ip;
        end;
        modelF(patchiter).inverse=Dout.inv{val}.inverse;
    end; % for patchiter
else %% RUN IN SERIAL - the code below should be the same as the code above -need to do this more elegantly
    for patchiter=1:Npatchiter,
        Din=Dtest{1};
        Dout=Din;
        Ip=allIp(patchiter).Ip;
        
        switch funcname,
            case 'Classic',
                Din.inv{val}.inverse.Ip=Ip;
                Dout    = spm_eeg_invert_classic(Din);
            case 'Current',
                warning('Patch centres are currently fixed for this algorithm (iteration will have no effect!)');
                Dout    = spm_eeg_invert(Din); %
                Dout.inv{val}.inverse.Ip=Ip;
        end;
        modelF(patchiter).inverse=Dout.inv{val}.inverse;
    end; % for patchiter
end; %% if par



for patchiter=1:Npatchiter,
    allF(patchiter)=modelF(patchiter).inverse.F;
    manyinverse{patchiter}=modelF(patchiter).inverse;
end;

[bestF,bestind]=max(allF);
disp('model evidences relative to maximum:')

sort(allF-bestF)

Dtest{1}.inv{val}.inverse=modelF(bestind).inverse; %% return best model for now
if Npatchiter>1, %% keep iterations if more than 1
%% commented out section for storing all iterations in dataset
    %     for f=1:Npatchiter,
%         Dtest{1}.inv{f+val}.inverse=modelF(f).inverse; %% set fields in inversion to specific iterations
%         Dtest{1}.inv{f+val}.comment=sprintf('Iteration %d of %d',f,Npatchiter);
%     end;
%     
    if (Dtest{1}.inv{val}.inverse.BMAflag==1)
        disp('Running BMA to get current estimate');
        [Jbma,qCbma]=spm_eeg_invert_bma(manyinverse,allF); %% onlt the mean is calculated using BMA (but could extend to covariance)
        %Dtest{1}.inv{val}.inverse.T=1; %% Jbma is the sum of all modes
        Dtest{1}.inv{val}.inverse.J={Jbma};
        Dtest{1}.inv{val}.inverse.qC=qCbma;
        Dtest{1}.inv{val}.inverse.allF=allF;
        %Dtest{1}.inv{val}.comment={sprintf('BMA of %d solutions',Npatchiter)};
    else % NOT BMA- just take the best
        disp('Using best patch set to current estimate');
        
        Dtest{1}.inv{val}.comment{1}=sprintf('Best F of %d solutions',Npatchiter);
       % keyboard
      %  [Dtest{1}.inv{1}.comment{1}]
        
        Dtest{1}.inv{val}.inverse=modelF(bestind).inverse; %% return best model for now
        Dtest{1}.inv{val}.inverse.allF=allF;
    end; % if BMA
end;








spm_eeg_invert_display(Dtest{1});









