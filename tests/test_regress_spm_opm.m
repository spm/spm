function tests = test_regress_spm_opm
% regresion test for OPM functions
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_regress_spm_opm_1(testCase)
spm('defaults','eeg');

data = fullfile(spm('Dir'),'tests','data','OPM','OPM_meg_001.cMEG');
positions = fullfile(spm('Dir'),'tests','data','OPM','OPM_HelmConfig.tsv');
%- read data 
%--------------------------------------------------------------------------
S = [];
S.data =data;
S.positions= positions;
D = spm_opm_create(S);


%- highpass
%--------------------------------------------------------------------------
S=[];
S.D=D;
S.freq=[10];
S.band = 'high';
fD = spm_eeg_ffilter(S);
delete(D);

%- lowpass
%--------------------------------------------------------------------------
S = [];
S.D = fD;
S.freq = [70];
S.band = 'low';
ffD = spm_eeg_ffilter(S);
delete(fD);

%- adaptive multipole models
%--------------------------------------------------------------------------
S = [];
S.D = ffD;
S.corrLim = .98;
mD = spm_opm_amm(S);
delete(ffD)

%- epoch 
%--------------------------------------------------------------------------
S =[];
S.D=mD;
S.timewin=[-50 200];
S.triggerChannels ={'Trigger 6 Z'};
eD= spm_opm_epoch_trigger(S);
delete(mD)

S=[];
S.D = eD;
S.timewin = [-50 0];
bD = spm_eeg_bc(S);
delete(eD);

%- average
%--------------------------------------------------------------------------
S=[];
S.D=bD;
muD = spm_eeg_average(S);
delete(bD);

MEGind = indchantype(eD,'MEGMAG');
used = setdiff(MEGind,badchannels(muD));
pl =muD(used,:,:)';

delete(muD);
peakValInRange = max(abs(pl(:))) < 154.6222 &  max(abs(pl(:))) > 154.3132;

testCase.verifyTrue(peakValInRange);