function [spatialmodename,Nmodes,newpctest,testchans]=spm_eeg_inv_prep_modes_xval(filenames, Nmodes, spatialmodename,Nblocks,pctest);
%function [spatialmodename,Nmodes,newpctest,testchans]=spm_eeg_inv_prep_modes_xval(filenames, Nmodes, spatialmodename,Nblocks,pctest);


%% prepare a spatial mode file for inversion
%% this file ensures the same spatial modes are used across different files (which would contain the same data but different head-models for example)
%% it also makes sure that the same channels groups are preserved to allow comparable cross validation and free energy metrics
%% input a list of the M/EEG dataset names: filenames
% Nmodes - number of required spatial modes (if empty uses all available
% channels)
% channels)
% spatialmodename- name of output file
% Nblocks- number of cross validation runs (optional and
% default 1)
% pctest- percentatge of channels to be used for testdata (optional and
% default 0)
%% if pctest*Nblocks=100 then will use independent MEG channels and may adjust pctest (in output) to 
%% accomodate this. ( k-fold cross val)
%% if pctest*Nblocks~=100 then uses random selection of pctest channels for each block (Monte-Carlo cross val)

%% in output file
%% megind- good meg channel indices
%% testchans - indices to megind of channels to be turned off during training phase (and tested later)
%% U{} - a different spatial modes matrix for each set of training channels or megind without indexed testchans or megind(setdiff(1:length(megind),testchans(b,:)))

% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
%
% Gareth Barnes
% $Id: spm_eeg_inv_prep_modes_xval.m 1 2016-09-12 15:02:25Z gareth $



if nargin<4,
    Nblocks=[];
end;
if nargin<5,
    pctest=[];
end;

if isempty(Nblocks),
    Nblocks=1;
end;
if isempty(pctest),
    pctest=0;
end;

Nfiles=size(filenames,1);
%% get lead fields for all files /headmodels and get unbiased mixture

D=spm_eeg_load(deblank(filenames(1,:)));
val=D.val;
megind=D.indchantype(D.inv{val}.forward.modality);
origbadchans=D.badchannels;

megind=setdiff(megind,origbadchans);
fprintf('Removed %d bad channels\n',length(origbadchans));

newpctest=pctest;
if round(Nblocks*pctest)==100,
    newpctest=floor(length(megind)/Nblocks)/length(megind)*100;
    fprintf('\nAdjusting pc test from %3.2f to %3.2f percent to make use of most MEG channels',pctest,newpctest);
    
else
    fprintf('\nTaking random selections of %3.2f percent of channels',newpctest);
end;
    

Ntest=round(newpctest*length(megind)/100); %% number of test channels per block

testchans=zeros(Nblocks,Ntest);
Ntrain=length(megind)-Ntest;

if isempty(Nmodes),
    Nmodes=length(megind)-Ntest;
    fprintf('\nUsing maximum of %d spatial modes\n',Nmodes);
else
    fprintf('\nUsing %d spatial modes\n',Nmodes);
end;
allL=spm_eeg_lgainmat(D);

if size(allL,1)~=length(megind),
    error('Mismatch in channel numbers (internal error)');
end;


U={};
useind=zeros(Nblocks,Ntrain);
testchans=zeros(Nblocks,Ntest);

alldum=randperm(length(megind));

for b=1:Nblocks,
    if round(Nblocks*pctest)==100,
        dum=alldum((b-1)*Ntest+1:b*Ntest);
    else
        dum=randperm(length(megind));
    end;

    testchans(b,:)=dum(1:Ntest); %% channels which will be used for testing (Set to bad during training)
    useind(b,:)=setdiff(1:length(megind),testchans(b,:)); %% training channels used here to get the spatial modes
end;

for f=2:Nfiles,
    D=spm_eeg_load(deblank(filenames(f,:)));
    val=D.val;
    if D.indchantype(D.inv{val}.forward.modality)~=megind,
        error('different active channels in these files');
    end;
    L=spm_eeg_lgainmat(D);
    allL=[allL L]; %% make composite lead field matrix
end; % for f


U={};

for b=1:Nblocks,
    fprintf('\nPreparing modes file %d block %d of %d for %d training and %d test chans\n',f, b,Nblocks,Ntrain,Ntest);
    if Nmodes<Ntrain,
        [U1,S,V]=spm_svd(allL(useind(b,:),:)*allL(useind(b,:),:)',1e-12); %% get general spatial modes matrix
        U{b}=U1(:,1:Nmodes)';
    else
        U{b}=eye(Ntrain);
    end; % if Nmodes
end; %for b

fprintf('\n saving spatial mode file %s\n',spatialmodename);
save(spatialmodename,'U','testchans','megind'); %%
