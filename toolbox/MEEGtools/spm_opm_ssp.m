function [sspD] = spm_opm_ssp(S)
% Removes eignevectors of training data covaraince from test data
% FORMAT D = spm_opm_ssp(S)
%   S               - input structure
%  fields of S:
%   S.trainD    - SPM MEEG object         - Default: no Default
%   S.testD     - SPM MEEG object         - Default: training data
%   S.trainwin  - Training Epoch (s)      - Default: [0,Inf]
%   S.testwin   - Testing Epoch (s)       - Default: [0,Inf]
%   S.ncomp     - Number of eigenvectors  - Default: 8
%   S.weights   - Channel weighting       - Default: ones(1,nchannels)
%   S.hp        - high pass filter (Hz)   - Default: no filter
%   S.lp        - low pass filter  (Hz)   - Default: no filter
%   S.notch     - Notch filter (Hz)       - Default: no filter
% Output:
%   D               - denoised MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright  Tim Tierney
%
% Example 
% S=[];
% S.trainD =D;                        % D object
% S.testD = D;                        % D object(potentially different);
% S.trainwin = [1,60];                % epochs 
% S.testwin = [61,Inf];               % epochs      
% S.ncomp = 32;                       % number of eigenvectors
% S.weights  =ones(1,length(megind)); % vector of weights(1 for each channel)
% S.weights(t2inds)= 1/6^2;           % e.g factor of six difference in noise floor
% S.lp=[100];                         % filter only pplied to training data
% S.hp=[1];                           % filter only pplied to training data
% S.notch = [47 53;97 103];           % filter only pplied to training data
% sspD = spm_opm_ssp(S);

%- Arg Check
%--------------------------------------------------------------------------
errorMsg = 'an training dataset must be supplied.';
if ~isfield(S, 'trainD'),             error(errorMsg); end
if ~isfield(S, 'testD'),              S.testD=S.trainD; end
if ~isfield(S, 'trainwin'),           S.trainwin=[0,Inf]; end
if ~isfield(S, 'testwin'),            S.testwin=[0,Inf]; end
if ~isfield(S, 'ncomp'),              S.ncomp=8; end
if ~isfield(S, 'hp'),                 S.hp=0; end
if ~isfield(S, 'lp'),                 S.lp=0; end
if ~isfield(S, 'notch'),              S.notch=0; end

%-Get datasets
%--------------------------------------------------------------------------

if isfield(S, 'trainwin')
  args =[];
  args.D= S.trainD;
  args.timewin = S.trainwin*1000;
  args.prefix = 'train_';
  Dtrain = spm_eeg_crop(args);
  Dtrain.save();
else
  Dtrain = S.train;
end

if isfield(S, 'testwin')
  args =[];
  args.D= S.testD;
  args.timewin = S.testwin*1000;
  args.prefix = 'ssp_';
  sspD = spm_eeg_crop(args);
  sspD.save();
else
  sspD = S.test;
end

%-Get usable channels
%--------------------------------------------------------------------------
chaninds = indchantype(Dtrain,'MEG');
badinds = badchannels(Dtrain);
usedinds = setdiff(chaninds,badinds);
usedLabs= chanlabels(Dtrain,usedinds);

%- dc Correct training data
%--------------------------------------------------------------------------

args=[];
args.D= Dtrain;
args.timewin = [Dtrain.time(1)*1000 Dtrain.time(end)*1000];
args.save = false;
Dtrain = spm_eeg_bc(args);
Dtrain.save();


%- filter training data
%--------------------------------------------------------------------------
if(S.hp>0)
  args =[];
  args.D=Dtrain;
  args.freq=S.hp;
  args.band = 'high';
  args.order = 3;
  args.save = false;
  Dtrain = spm_eeg_ffilter(args);
end

if(S.lp>0)
  args = [];
  args.D = Dtrain;
  args.freq = S.lp;
  args.band = 'low';
  args.order = 6;
  args.save = false;
  Dtrain = spm_eeg_ffilter(args);
end

if(all(S.notch>0))
  for i = 1:size(S.notch,1)
    args = [];
    args.D = Dtrain;
    args.freq = S.notch(i,:);
    args.band = 'stop';
    args.order = 2;
    args.save = false;
    Dtrain = spm_eeg_ffilter(args);
  end
end

%- Define Projection
%--------------------------------------------------------------------------
chunkSize = round(Dtrain.fsample);
C2 = zeros(length(usedinds),length(usedinds));
N= size(Dtrain,2);

for startIdx = (3*chunkSize):chunkSize:(N-3*chunkSize)
  endIdx = min(startIdx + chunkSize - 1, N-3*chunkSize);
  tmp = Dtrain(:,startIdx:endIdx, :)';
  chunk = tmp(:,usedinds);
  C2 = C2 + (chunk'*chunk);
end

Dtrain.delete();

[U,~,~] = svd(C2);
Us = U(:,1:S.ncomp);
W = eye(size(Us,1));

for i = 1:size(Us,1)
    W(i,i)=S.weights(i);
end

Ui = pinv(Us'*W*Us)*Us'*W;
M = eye(size(Us,1))-Us*Ui;


%- Work out chunk size
%--------------------------------------------------------------------------
begs=1:(chunkSize*10):size(sspD,2);
ends = (begs+chunkSize*10-1);
if(ends(end)>size(sspD,2))
    ends(end)= size(sspD,2);
end

%- Apply Projeciton
%--------------------------------------------------------------------------

for i =1:length(begs)
  inds = begs(i):ends(i);
  % sspD(Yinds,inds,j)=M*S.D(Yinds,inds,j) is slow (disk read)
  out = sspD(:,inds,:);
  Y=out(usedinds,:);
  out(usedinds,:)=M*Y;
  sspD(:,inds,:)=out;
end
sspD.save();

%-Update forward modelling information
%--------------------------------------------------------------------------
grad = sspD.sensors('MEG');
[~,sinds] = spm_match_str(usedLabs,grad.label);
tmpTra= eye(size(grad.coilori,1));
tmpTra(sinds,sinds)=M;
grad.tra= tmpTra*grad.tra;
grad.balance.previous= grad.balance.current;
grad.balance.current    = 'ssp';
sspD = sensors(sspD,'MEG',grad);
if isfield(sspD,'inv')
    if isfield(sspD.inv{1},'gainmat')
        fprintf(['Clearing current forward model, please recalculate '...
            'with spm_eeg_lgainmat\n']);
        sspD.inv{1} = rmfield(sspD.inv{1},'gainmat');
    end
    if isfield(sspD.inv{1},'datareg')
        sspD.inv{1}.datareg.sensors = grad;
    end
    if isfield(sspD.inv{1},'forward')
        voltype = sspD.inv{1}.forward.voltype;
        sspD.inv{1}.forward = [];
        sspD.inv{1}.forward.voltype = voltype;
        sspD = spm_eeg_inv_forward(sspD,1);
    end
end
sspD.save();
end
