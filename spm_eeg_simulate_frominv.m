function [Dnew]=spm_eeg_simulate_frominv(D,prefix,val,whitenoise,SNRdB,trialind)
% Project a source inversion solution back out to the sensor level plus some noise
% D          - original dataset
% prefix     - prefix of new dataset
% val        - use solution (and lead fields) corresponding to this index
% whitenoise - total rms white noise in Tesla
% SNRdB      - SNR in dBs (alternative to specifying white noise)
% trialind   - trials in which the simulated signal is to appear
%              (all other trials will be noise)
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_eeg_simulate_frominv.m 5664 2013-10-01 18:39:05Z spm $


%-Load in original data
%==========================================================================


if nargin<2, prefix     = ''; end
if nargin<3, val        = []; end
if nargin<4, whitenoise = []; end
if nargin<5, SNRdB      = []; end
if nargin<6, trialind   = []; end

if isempty(prefix), prefix='sim'; end

if isempty(val)
    val=D{useind}.val;
end

useind = 1; % D to use

if ~xor(isempty(whitenoise),isempty(SNRdB))
    error('Must specify either white noise level or sensor level SNR');
end

newfilename = spm_file(D{useind}.fname,'prefix',prefix);
try
    Dnew = spm_eeg_load(newfilename);
    disp('Overwriting data in existing file');
catch
    Dnew = D{useind}.clone(newfilename);
end

if isempty(trialind)
    trialind = 1:Dnew.ntrials;
end

disp('Simulating data on MEG channels only'); % for now
chanind = strmatch('MEG',Dnew.chantype);


%-Some default noise levels
%==========================================================================
% 
% if isempty(whitenoise)
%     sensor_noise_TrtHz=10e-15; %% Sensor noise in Tesla per root Hz; default 10 fT/rtHz
%     sensor_bw_Hz=80; %% recording bandwith in Hz
%     whitenoise=sqrt(sensor_bw_Hz)*sensor_noise_TrtHz;
%     disp('setting default 10ftrtHz white noise in 80Hz BW');
% else
%     whitenoise=1;
% end


%-Create a new forward model
% (not necessarily the same as the one used to make the inversion)
%==========================================================================
fprintf('Computing Gain Matrix: ')
spm_input('Creating gain matrix',1,'d');    % Shows gain matrix computation

M     = Dnew.inv{val}.inverse.M;
J     = Dnew.inv{val}.inverse.J{1};
%Y    = Dnew.inv{val}.inverse.Y;
Lorig = Dnew.inv{val}.inverse.L;
T     = Dnew.inv{val}.inverse.T;
U     = Dnew.inv{val}.inverse.U;
%Is   = Dnew.inv{val}.inverse.Is;
Ic    = Dnew.inv{val}.inverse.Ic;
It    = Dnew.inv{val}.inverse.It;
try
    [L,Dnew] = parfor_spm_eeg_lgainmat(Dnew);  % Gain matrix- from file rather than from inversion itself
catch
    [L,Dnew] = spm_eeg_lgainmat(Dnew);         % Gain matrix- from file rather than from inversion itself
end
        
Lnew  = U*L;
if max(max(Lnew-Lorig))>0
    warning('Lead fields have changed since inversion (but keeping original U)');
end
        
tmp   = (Lnew*J)*T';
tmp   = U'*tmp; %% channels x time
% 
%tmpo = (Lorig*J)*T';
%tmpo = U'*tmpo;

switch Dnew.sensors('MEG').chanunit{1}
    case 'T'
        whitenoise = whitenoise; %% rms tesla

    case 'fT'
        whitenoise = whitenoise*1e15; %% rms femto tesla

    otherwise
        error('unknown sensor unit')
end


allchanstd    = std(tmp');
meanrmssignal = mean(allchanstd);


if ~isempty(SNRdB)
    whitenoise = meanrmssignal.*(10^(SNRdB/20));
    disp(sprintf('Setting white noise to give sensor level SNR of %dB',SNRdB));
end


for i=1:Dnew.ntrials
    Dnew(:,:,i) = Dnew(:,:,i).*0;
    if find(i,trialind), %% only add signal to specific trials
        Dnew(Ic,It,i) = tmp;
    else
        Dnew(Ic,It,i) = zeros(size(tmp));
    end
    Dnew(:,:,i) = Dnew(:,:,i) + randn(size(Dnew(:,:,i))).*whitenoise; %% add white noise in fT
end


%-Plot and save
%==========================================================================
[dum,plotind] = sort(allchanstd);

% vert=Dnew.inv{val}.mesh.tess_mni.vert;
% 
% Nj        = size(vert,1);
% M     = mean(X(:,f1ind)'.^2,1);
% G     = sqrt(sparse(1:Nj,1,M,Nj,1));
% Fgraph    = spm_figure('GetWin','Graphics');
% j     = find(G);
% 
% clf(Fgraph)
% figure(Fgraph)
% spm_mip(G(j),vert(j,:)',6);
% axis image
% title({sprintf('Generated source activity')});
% drawnow

figure
hold on
aux = tmp(plotind(end),:);
subplot(2,1,1);
plot(Dnew.time,Dnew(plotind(end),:,trialind(1)),'b',Dnew.time(It),aux,'r');
title('Measured activity over max sensor');
legend('Noisy','Noiseless');
subplot(2,1,2);
aux = tmp(plotind(floor(length(plotind)/2)),:);
plot(Dnew.time,Dnew(plotind(floor(length(plotind)/2)),:,trialind(1)),'b',Dnew.time(It),aux,'r');
title('Measured activity over median sensor');
legend('Noisy','Noiseless');

Dnew.save;

fprintf('\n Finish\n')
