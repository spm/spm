function [Dnew,meshsourceind]=spm_eeg_simulate(D,prefix,patchmni,simsignal,ormni,woi,whitenoise,SNRdB,trialind,mnimesh,SmthInit);
%function [Dnew,meshsourceind]=spm_eeg_simulate(D,prefix,patchmni,simsignal,woi,whitenoise,SNRdB,trialind,mnimesh,SmthInit);
%% Simulate a number of MSP patches at specified locations on existing mesh
%
% Created by:   Jose David Lopez - ralph82co@gmail.com
%               Gareth Barnes - g.barnes@ucl.ac.uk
%               Vladimir Litvak - litvak.vladimir@gmail.com
%
%% D dataset
%% prefix : prefix of new simulated dataset
%% patchmni : patch centres in mni space or patch indices
%% simsignal : Nsourcesx time series withinn woi
%% woi: window of interest in seconds
%% whitenoise level in Tesla
%% SNRdB power signal to noise ratio in dBs
%% trialind: trials on which the simulated data will be added to the noise
%% mnimesh : a new mesh with vertices in mni space
%% SmthInit - the smoothing step that creates the patch- larger numbers larger patches default 0.6. Note current density should be constant (i.e. larger patch on tangential surface will not give larger signal)
%% Outputs
%% Dnew- new dataset
%% meshsourceind- vertex indices of sources on the mesh
% $Id: spm_eeg_simulate.m 6077 2014-06-30 16:55:03Z spm $

%% LOAD IN ORGINAL DATA
useind=1; % D to use
if nargin<2,
    prefix='';
end;

if nargin<3,
    patchmni=[];
end;


if nargin<4,
    simsignal=[];
end;
if nargin<5,
    ormni=[];
end;

if nargin<6,
    woi=[];
end;


if nargin<7,
    whitenoise=[];
end;

if nargin<8,
    SNRdB=[];
end;

if nargin<9,
    trialind=[];
end;

if nargin<10,
    mnimesh=[];
end;



if nargin<11
    SmthInit=[]; %% number of iterations used to smooth patch out (more iterations, larger patch)
end;



if isempty(prefix),
    prefix='sim';
end;

if isempty(SmthInit),
    SmthInit=0.6;
end;

if isempty(woi),
    woi=[D{useind}.time(1) D{useind}.time(end)];
end;


val=D{useind}.val;



if isempty(patchmni),
    patchmni=[-45.4989  -30.6967    4.9213;...
        46.7322  -31.2311    4.0085];
    
end;


if ~xor(isempty(whitenoise),isempty(SNRdB))
    error('Must specify either white noise level or sensor level SNR');
end;





[a1 b1 c1]=fileparts(D{useind}.fname);
newfilename=[prefix b1];

%% forcing overwrite of an existing file
Dnew=D{useind}.clone([prefix b1]);


if isempty(trialind)
    trialind=1:Dnew.ntrials;
end;

modstr=deblank(modality(D{1}));
disp(sprintf('Simulating data on %s channels only',modstr));



if ~isempty(mnimesh),
    Dnew.inv{val}.mesh.tess_mni.vert=mnimesh.vert;
    Dnew.inv{val}.mesh.tess_mni.face=mnimesh.face;
    Dnew.inv{val}.forward.mesh.vert=spm_eeg_inv_transform_points(Dnew.inv{val}.datareg.fromMNI,mnimesh.vert);
    Dnew.inv{val}.forward.mesh.face=mnimesh.face;
end; % if

% Two synchronous sources
if patchmni~=0,
    Ndips=size(patchmni,1);
else
    Ndips=0;
end;

if size(simsignal,1)~=Ndips,
    error('number of signals given does not match number of sources');
end;

meshsourceind=[];

disp('Using closest mesh vertices to the specified coordinates')
for d=1:Ndips,
    vdist= Dnew.inv{val}.mesh.tess_mni.vert-repmat(patchmni(d,:),size(Dnew.inv{val}.mesh.tess_mni.vert,1),1);
    dist=sqrt(dot(vdist',vdist'));
    [mnidist(d),meshsourceind(d)] =min(dist);
end;

disp(sprintf('Furthest distance %3.2f mm',max(mnidist)));
if max(mnidist)>0.1
    warning('Supplied vertices do not sit on the mesh!');
end;


Ndip = size(simsignal,1);       % Number of dipoles




try Dnew.inv{val}.forward.vol.unit
    switch(Dnew.inv{val}.forward.vol.unit), %% correct for non-SI lead field scaling
        case 'mm'
            
            Lscale=1000*1000;
        case 'cm'
            
            Lscale=100*100;
        case 'm'
            Lscale=1.0;
            
        otherwise
            error('unknown volume unit');
            
    end;
catch
    disp('No units found');
    Lscale=1.0;
end;


%% WAVEFORM FOR EACH SOURCE

Ntrials = Dnew.ntrials;             % Number of trials

% define period over which dipoles are active
startf1  = woi(1);                  % (sec) start time
endf1 = woi(2); %% end time
f1ind = intersect(find(Dnew.time>startf1),find(Dnew.time<=endf1));

if length(f1ind)~=size(simsignal,2),
    error('Signal does not fit in time window');
end;


if isequal(modstr, 'MEG')
    chanind = Dnew.indchantype({'MEG', 'MEGPLANAR'}, 'GOOD');
else
    chanind = Dnew.indchantype(modality, 'GOOD');
end

units = Dnew.units;
indunits=Dnew.indchannel(Dnew.inv{val}.forward.sensors.label);
units=units(indunits);

labels=Dnew.chanlabels(chanind);


chans = Dnew.indchantype(modstr, 'GOOD');



if ~isempty(ormni), %%%% DIPOLE SIMULATION
    disp('SIMULATING DIPOLE SOURCES');
    if size(ormni)~=size(patchmni),
        error('A 3D orientation must be specified for each source location');
    end;
    
    posdipmm=Dnew.inv{val}.datareg.fromMNI*[patchmni ones(size(ormni,1),1)]'; %% put into MEG space
    posdipmm=posdipmm(1:3)';
    %% need to make a pure rotation for orientation transform to native space
    M1 = Dnew.inv{val}.datareg.fromMNI;
    [U, L, V] = svd(M1(1:3, 1:3)); 
    ordip=ormni*(U*V');
    ordip=ordip./sqrt(dot(ordip,ordip)); %% make sure it is unit vector
    
    %% NB COULD ADD A PURE DIPOLE SIMULATION IN FUTURE
    sens=Dnew.inv{val}.forward.sensors;
    vol=Dnew.inv{val}.forward.vol;
    
    
    
    
    ftchans=Dnew.inv{val}.forward.sensors.label;
    for j=1:numel(chans),
        usedind(j)=strmatch(labels{j},ftchans);
    end;
    
    
    tmp=zeros(length(chanind),Dnew.nsamples);
    for i=1:Ndip,
        gmn = ft_compute_leadfield(posdipmm(i,:)*1e-3, sens, vol,  'dipoleunit', 'nA*m','chanunit',units);
        gain=gmn*ordip';
        tmp(:,f1ind)=tmp(:,f1ind)+gain(usedind,:)*simsignal(i,:);
    end; % for i
    
    
else %%% CURRENT DENSITY ON SURFACE SIMULATION
    disp('SIMULATING CURRENT DISTRIBUTIONS ON MESH');
    %% CREATE A NEW FORWARD model for e mesh
    fprintf('Computing Gain Matrix: ')
    spm_input('Creating gain matrix',1,'d');    % Shows gain matrix computation
    
    [L Dnew] = spm_eeg_lgainmat(Dnew);              % Gain matrix
    
    Nd    = size(L,2);                          % number of dipoles
    X     = zeros(Nd,size(Dnew,2));                     % Matrix of dipoles
    fprintf(' - done\n')
    
    
    
    % Green function for smoothing sources with the same distribution than SPM8
    fprintf('Computing Green function from graph Laplacian:')
    
    
    vert  = Dnew.inv{val}.mesh.tess_mni.vert;
    face  = Dnew.inv{val}.mesh.tess_mni.face;
    A     = spm_mesh_distmtx(struct('vertices',vert,'faces',face),0);
    GL    = A - spdiags(sum(A,2),0,Nd,Nd);
    GL    = GL*SmthInit/2;
    Qi    = speye(Nd,Nd);
    QG    = sparse(Nd,Nd);
    for i = 1:8,
        QG = QG + Qi;
        Qi = Qi*GL/i;
    end
    QG    = QG.*(QG > exp(-8));
    QG    = QG*QG;
    %QG=QG./sum(sum(QG));
    clear Qi A GL
    fprintf(' - done\n')
    
    
    % Add waveform of all smoothed sources to their equivalent dipoles
    % QGs add up to 0.9854
    fullsignal=zeros(Ndip,Dnew.nsamples); %% simulation padded with zeros
    fullsignal(1:Ndip,f1ind)=simsignal;
    for j=1:Ndip
        for i=1:Dnew.nsamples,
            X(:,i) = X(:,i) + fullsignal(j,i)*QG(:,meshsourceind(j)); %% this will be in Am
        end
    end
   
    % Copy same data to all trials
    tmp=L*X;
    
end; % if ori

if isfield(Dnew.inv{val}.forward,'scale'),
        tmp=Lscale.*tmp./Dnew.inv{val}.forward.scale; %% account for rescaling of lead fields
    else
        tmp=Lscale.*tmp; %% no rescaling
    end;
    
    
try
    
    switch Dnew.sensors(modstr).chanunit{1}
        case 'T'
            whitenoise=whitenoise; %% rms tesla
            tmp=tmp;
        case 'fT'
            whitenoise=whitenoise*1e15; %% rms femto tesla
            tmp=tmp*1e15;
        case 'V'
            whitenoise=whitenoise; %% volts
            tmp=tmp;
        case 'uV'
            whitenoise=whitenoise*1e6; %% micro volts
            tmp=tmp;
            
            
        otherwise
            error('unknown sensor unit')
    end;
    
catch
    disp('No sensor units found');
end;



allchanstd=std(tmp');
meanrmssignal=mean(allchanstd);


if ~isempty(SNRdB),
    whitenoise = meanrmssignal.*(10^(-SNRdB/20));
    disp(sprintf('Setting white noise to give sensor level SNR of %dB',SNRdB));
end;



for i=1:Ntrials
    if any(i == trialind), %% only add signal to specific trials
        Dnew(chans,:,i) = tmp;
    else
        Dnew(chans,:,i)=zeros(size(tmp));
    end;
    Dnew(:,:,i)=Dnew(:,:,i)+randn(size(Dnew(:,:,i))).*whitenoise; %% add white noise in fT
end




%% Plot and save
[dum,tmpind]=sort(allchanstd);
dnewind=chans(tmpind);

if isempty(ormni)
    Nj      = size(vert,1);
    M       = mean(X(:,f1ind)'.^2,1);
    G       = sqrt(sparse(1:Nj,1,M,Nj,1));
    Fgraph  = spm_figure('GetWin','Graphics');
    j       = find(G);
    
    clf(Fgraph)
    figure(Fgraph)
    spm_mip(G(j),vert(j,:)',6);
    axis image
    title({sprintf('Generated source activity')});
    drawnow
end;
figure
hold on

aux = tmp(tmpind(end),:);
subplot(2,1,1);
plot(Dnew.time,Dnew(dnewind(end),:,1),Dnew.time,aux,'r');
title('Measured activity over max sensor');
legend('Noisy','Noiseless');
subplot(2,1,2);
aux = tmp(tmpind(floor(length(tmpind)/2)),:);
plot(Dnew.time,Dnew(dnewind(floor(length(tmpind)/2)),:,1),Dnew.time,aux,'r');
title('Measured activity over median sensor');
legend('Noisy','Noiseless');

Dnew.save;

fprintf('\n Finish\n')

