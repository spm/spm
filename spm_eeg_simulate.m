function [Dnew,meshsourceind,signal]=spm_eeg_simulate(D,prefix,patchmni,dipfreq,woi,dipmoment,whitenoise,SNRdB,trialind,mnimesh,SmthInit);
%% Simulate a number of MSP patches at specified locations on existing mesh
% Synthetic MEG data generator for SPM8
% This is a demo version related with the Technical Note:
% XXX
%
% Created by:	Jose David Lopez - ralph82co@gmail.com
%				Gareth Barnes - g.barnes@ucl.ac.uk
%				Vladimir Litvak - litvak.vladimir@gmail.com
%
%% D dataset
%% prefix : prefix of new simulated dataset
%% patchmni : patch centres in mni space or patch indices
%% dipfreq : frequency of simulated sources (Hz)
%% mnimesh : a new mesh with vertices in mni space
%% dipmoment: relative dipole moments
%% woi : time window of source activity
%% SmthInit - the smoothing step that creates the patch- larger numbers larger patches default 0.6. Note current density should be constant (i.e. larger patch on tangential surface will not give larger signal)
%
% $Id: spm_eeg_simulate.m 5478 2013-05-08 15:00:56Z gareth $

%% LOAD IN ORGINAL DATA
useind=1; % D to use
if nargin<2,
    prefix='';
end;

if nargin<3,
    patchmni=[];
end;


if nargin<4,
    freqs=[];
end;
if nargin<5,
    woi=[];
end;

if nargin<6,
    dipmoment=[];
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

if ~isempty(SNRdB),
    dipmoment=ones(size(patchmni,1),1).*20e-9; %% set to 20nAm
end;




[a1 b1 c1]=fileparts(D{useind}.fname);
newfilename=[prefix b1];
try
    Dnew=spm_eeg_load(newfilename);
    disp('Overwriting data in existing file');
catch
    Dnew=D{useind}.clone([prefix b1]);
end;

if isempty(trialind)
    trialind=1:Dnew.ntrials;
end;
disp('Simulating data on MEG channels only'); % for now
chanind=strmatch('MEG',Dnew.chantype);

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

if length(dipfreq)~=Ndips,
    error('number of frequencies given does not match number of sources');
end;

meshsourceind=[];

disp('Using closest mesh vertices to the specified coordinates')
for d=1:Ndips,
    vdist= Dnew.inv{val}.mesh.tess_mni.vert-repmat(patchmni(d,:),size(Dnew.inv{val}.mesh.tess_mni.vert,1),1);
    dist=sqrt(dot(vdist',vdist'));
    [mnidist(d),meshsourceind(d)] =min(dist);
end;

disp(sprintf('Furthest distance %3.2f',max(mnidist)));
if max(mnidist)>0.1
    warning('Supplied vertices do not sit on the mesh!');
end;


Ndip = size(meshsourceind,2);		% Number of dipoles
if isempty(dipfreq),
    dipfreq = ones(Ndips,1).*20;					% Source frequency
end;

%% some default noise levels
% 
% if isempty(whitenoise)&&isempty(SNRdB),
%     sensor_noise_TrtHz=10e-15; %% Sensor noise in Tesla per root Hz; default 10 fT/rtHz
%     sensor_bw_Hz=80; %% recording bandwith in Hz
%     whitenoise=sqrt(sensor_bw_Hz)*sensor_noise_TrtHz;
%     disp('setting default 10ftrtHz white noise in 80Hz BW');
% else
%     whitenoise=1;
% end;




switch(Dnew.inv{val}.forward.vol.unit), %% correct for non-SI lead field scaling
    case 'mm'
        
        Lscale=1000*1000;
    case 'cm'
        
        Lscale=100*100;
    otherwise
        error('unknown volume unit')
end;


%% WAVEFORM FOR EACH SOURCE

Ntrials = Dnew.ntrials;				% Number of trials

% define period over which dipoles are active
startf1  = woi(1);					% (sec) start time
endf1 = woi(2); %% end time
f1ind = intersect(find(Dnew.time>startf1),find(Dnew.time<=endf1));

% Create the waveform for each source
signal = zeros(Ndip,length(Dnew.time));
for j=1:Ndip				% For each source
    for i=1:Ntrials			% and trial
        f1 = dipfreq(j);	% Frequency depends on stim condition
        amp1 = dipmoment(j);	% also the amplitude
        phase1 = pi/2;
        signal(j,f1ind) = signal(j,f1ind)...
            + amp1*sin((Dnew.time(f1ind)...
            - Dnew.time(min(f1ind)))*f1*2*pi + phase1);
    end
end

%% CREATE A NEW FORWARD PROBLEM

fprintf('Computing Gain Matrix: ')
spm_input('Creating gain matrix',1,'d');	% Shows gain matrix computation

[L Dnew] = spm_eeg_lgainmat(Dnew);				% Gain matrix

Nd    = size(L,2);							% number of dipoles
X	  = zeros(Nd,size(Dnew,2));						% Matrix of dipoles
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
for j=1:Ndip
    for i=1:size(Dnew,2)
        X(:,i) = X(:,i) + signal(j,i)*QG(:,meshsourceind(j)); %% this will be in Am
    end
end


% Copy same data to all trials
tmp=L*X;
if isfield(Dnew.inv{val}.forward,'scale'),
    tmp=Lscale.*tmp./Dnew.inv{val}.forward.scale; %% account for rescaling of lead fields
else
    tmp=Lscale.*tmp; %% no rescaling
end; 

switch Dnew.sensors('MEG').chanunit{1}
    case 'T'
        whitenoise=whitenoise; %% rms tesla
        tmp=tmp;
    case 'fT'
        whitenoise=whitenoise*1e15; %% rms femto tesla
        tmp=tmp*1e15;
    otherwise
        error('unknown sensor unit')
end;


allchanstd=std(tmp');
meanrmssignal=mean(allchanstd);


if ~isempty(SNRdB),
    whitenoise = meanrmssignal.*(10^(SNRdB/20));
    disp(sprintf('Setting white noise to give sensor level SNR of %dB',SNRdB));
end;


for i=1:Ntrials
    if any(i == trialind), %% only add signal to specific trials
        Dnew(:,:,i) = tmp;
    else
        Dnew(:,:,i)=zeros(size(tmp));
    end;
    Dnew(:,:,i)=Dnew(:,:,i)+randn(size(Dnew(:,:,i))).*whitenoise; %% add white noise in fT
end




%% Plot and save
[dum,plotind]=sort(allchanstd);


Nj		= size(vert,1);
M		= mean(X(:,f1ind)'.^2,1);
G		= sqrt(sparse(1:Nj,1,M,Nj,1));
Fgraph	= spm_figure('GetWin','Graphics');
j		= find(G);

clf(Fgraph)
figure(Fgraph)
spm_mip(G(j),vert(j,:)',6);
axis image
title({sprintf('Generated source activity')});
drawnow

figure
hold on
aux = tmp(plotind(end),:);
subplot(2,1,1);
plot(Dnew.time,Dnew(plotind(end),:,1),Dnew.time,aux,'r');
title('Measured activity over max sensor');
legend('Noisy','Noiseless');
subplot(2,1,2);
aux = tmp(plotind(floor(length(plotind)/2)),:);
plot(Dnew.time,Dnew(plotind(floor(length(plotind)/2)),:,1),Dnew.time,aux,'r');
title('Measured activity over median sensor');
legend('Noisy','Noiseless');

Dnew.save;

fprintf('\n Finish\n')

