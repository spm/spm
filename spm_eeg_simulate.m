function [D,meshsourceind,signal]=spm_eeg_simulate(D,prefix,patchmni,dipfreq,woi,dipmoment,mnimesh,keep_leads,SmthInit);
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
%% keep_leads: a flag if set to 0 then lead fields are re-calculated
%% dipmoment : dipole moment in Am
%% woi : time window of source activity
%% SmthInit - the smoothing step that creates the patch- larger numbers larger patches default 0.6. Note current density should be constant (i.e. larger patch on tangential surface will not give larger signal)
% 
% $Id: spm_eeg_simulate.m 5229 2013-02-04 14:02:29Z gareth $

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
    mnimesh=[];
end;

if nargin<8, %% use previously computed lead fields
    keep_leads=[];
end;

if nargin<9
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

if isempty(keep_leads),
    keep_leads=0;
end;
val=D{useind}.val;

if ~keep_leads,
    D{useind}.inv{val}.gainmat=''; %% force re-calculation of leads
end;

if isempty(patchmni),
         patchmni=[-45.4989  -30.6967    4.9213;...
   46.7322  -31.2311    4.0085];
    
end;

if isempty(dipmoment),
    dipmoment= ones(Ndips,1)*20*1e-9;		% dipole moment in Am
end;


[a1 b1 c1]=fileparts(D{useind}.fname);

Dnew=D{useind}.clone([prefix b1]);

disp('Simulating data on MEG channels only'); % for now
chanind=strmatch('MEG',Dnew.chantype);

if ~isempty(mnimesh),
    D{useind}.inv{1}.mesh.tess_mni.vert=mnimesh.vert;
    D{useind}.inv{1}.mesh.tess_mni.face=mnimesh.face;
    D{useind}.inv{1}.forward.mesh.vert=spm_eeg_inv_transform_points(D.inv{D.val}.datareg.fromMNI,mnimesh.vert);
    D{useind}.inv{1}.forward.mesh.face=mnimesh.face;
    
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

for d=1:Ndips,
    vdist= Dnew.inv{val}.mesh.tess_mni.vert-repmat(patchmni(d,:),size(Dnew.inv{val}.mesh.tess_mni.vert,1),1);
    dist=sqrt(dot(vdist',vdist'));
    [mnidist,meshsourceind(d)] =min(dist);
end;

Ndip = size(meshsourceind,2);		% Number of dipoles
if isempty(dipfreq),
    dipfreq = ones(Ndips,1).*20;					% Source frequency
end;

%% some default noise levels

sensor_noise_TrtHz=10e-15; %% Sensor noise in Tesla per root Hz; default 10 fT/rtHz
sensor_bw_Hz=80; %% recording bandwith in Hz
whitenoise=sqrt(sensor_bw_Hz)*sensor_noise_TrtHz;




 switch(Dnew.inv{val}.forward.vol.unit),
         case 'mm'
             dipmoment=dipmoment.*1000; %% Ampere mm
         case 'cm'
             dipmoment=dipmoment.*100; % Ampere cm
         otherwise
         error('unknown volume unit')
 end;

switch Dnew.sensors('MEG').chanunit{1}
    case 'T'
        whitenoise=whitenoise; %% rms tesla 
    case 'fT'
        whitenoise=whitenoise*1e15; %% rms femto tesla
     otherwise
        error('unknown sensor unit')
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
		X(:,i) = X(:,i) + signal(j,i)*QG(:,meshsourceind(j));
	end
end


% Copy same data to all trials 
tmp=L*X; 
tmp=tmp./repmat(Dnew.inv{val}.forward.scale,1,size(tmp,2)); %% account for scaling of lead fields

%meanrmssignal=mean(allchanstd);
for i=1:Ntrials
	Dnew(:,:,i) = tmp;
    Dnew(:,:,i)=Dnew(:,:,i)+randn(size(Dnew(:,:,i))).*whitenoise; %% add white noise in fT
end



%% Plot and save

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
aux = tmp(1,:);
plot(Dnew.time,Dnew(1,:,1),Dnew.time,aux,'r');
title('Measured activity over first sensor');
legend('Noisy','Noiseless');

Dnew.save;

fprintf('\n Finish\n')

