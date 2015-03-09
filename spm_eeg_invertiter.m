function [Dtest,modelF,allF]=spm_eeg_invertiter(Dtest,Npatchiter,funcname,patchind)

%  Function to perform several MSP type inversions with different
%  pseudo-randomly selected priors- in this case single cortical patches
%
% Npatchiter: number of iterations
% funcname is name of MSP alogorithm: current (spm_eeg_invert) or classic (spm_eeg_invert_classic)
% patchind is an optional list of indices of vertices which will be patch
% centres. patchind will have size Npatchiter*Np (where Np is number of patches set in
% inverse.Np )
% if  Dtest{1}.inv{val}.inverse.mergeflag==1 then merges Npatchiter posterior current
% distributions, else replaces posterior with best of the iterations.
% __________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
%
% Gareth Barnes
% $Id: spm_eeg_invertiter.m 6367 2015-03-09 15:02:25Z gareth $

if nargin<2,
    Npatchiter=[];
end;


if nargin<3,
    error('No function call specified');
end;

if nargin<4,
    patchind=[];
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





fprintf('Checking leadfields');


[L Dtest{1}] = spm_eeg_lgainmat(Dtest{1});  % Generate/load lead field- this stops it being done at each iteration


if isempty(patchind),
    disp('Reseting random number seed ! and then generating random patch centres ');
    rand('state',0);
    
    %% make random patch centers, except on the first iteration when we keep to a fixed stepsize
    for f=1:Npatchiter,
        tmp=randperm(Nvert);
        allIp(f).Ip=tmp(1:Np);
    end;
    allIp(1).Ip=[]; %% fix the step size of the first patch set (uses default MSP indices)
else
    for f=1:Npatchiter,
        allIp(f).Ip=patchind(f,:);
    end;
end;


modelF=[];
bestF=-Inf;
for patchiter=1:Npatchiter,
    %Din=Dtest{1};
    %Dout=Din;
    Ip=allIp(patchiter).Ip;
    disp(sprintf('patchiter %d/%d',patchiter,Npatchiter));
    switch funcname,
        case 'Classic',
            Dtest{1}.inv{val}.inverse.Ip=Ip;
            Dtest{1}    = spm_eeg_invert_classic(Dtest{1});
        case 'Current',
            warning('Patch centres are currently fixed for this algorithm (iteration will have no effect!)');
            Dtest{1}   = spm_eeg_invert(Dtest{1}); %
            %Dout.inv{val}.inverse.Ip=Ip;
        otherwise
            error('unknown function');
    end;
    if (Dtest{1}.inv{val}.inverse.mergeflag==1), %% will need all posteriors to merge them (could consider getting rid of worst ones here to save memory
        modelF(patchiter).inverse=Dtest{1}.inv{val}.inverse;
    else %% just keep track of best inversion
        if Dtest{1}.inv{val}.inverse.F>bestF,
            modelF(1).inverse=Dtest{1}.inv{val}.inverse;
            bestF=Dtest{1}.inv{val}.inverse.F;
        end;
    end;
    allF(patchiter)=Dtest{1}.inv{val}.inverse.F;
end; % for patchiter


%% will use these to merge posteriors or take best one



[bestF,bestind]=max(allF);
disp('model evidences relative to maximum:')

sort(allF-bestF)


%% for 1 iteration or for best solution
Dtest{1}.inv{val}.inverse=modelF(1).inverse; %% return best model (only 1 model saved in this option)
Dtest{1}.inv{val}.inverse.allF=allF;

if Npatchiter>1, %% keep iterations if more than 1
    
    if (Dtest{1}.inv{val}.inverse.mergeflag==1)
        Dtest{1}.inv{val}.inverse=modelF(bestind).inverse; % start with best model so far
        Qpriors=sparse(zeros(Npatchiter,size(modelF(1).inverse.qC,1)));
        for patchiter=1:Npatchiter,
            Qpriors(patchiter,:)=modelF(patchiter).inverse.qC;
        end;
        disp('Merging posterior distributions..');
        ugainfiles=Dtest{1}.inv{val}.gainmat;
        surfind=ones(Npatchiter,1);
        %% now mix the posteriors
        [Dtest{1}] = spm_eeg_invert_classic_mix(Dtest{1},val,Qpriors,surfind,ugainfiles);
        
        
        Dtest{1}.inv{val}.comment{1}=sprintf('Merged posterior of %d solutions',Npatchiter);
        mixF=Dtest{1}.inv{val}.inverse.F;
        if mixF-bestF<0
            warning('Merged solution is worse than the best');
        end;
        disp(sprintf('Improvement (when +ve) over best iteration %3.2f',mixF-bestF));
        
    else % NO merge - just take the best
        
        disp('Using best patch set to current estimate');
        Dtest{1}.inv{val}.comment{1}=sprintf('Best F of %d solutions',Npatchiter);
    end; % if BMA
    
    
end;








spm_eeg_invert_display(Dtest{1});









