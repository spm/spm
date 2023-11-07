function tests = test_spm_eeg_ffilter
% Unit Tests for spm_eeg_ffilter
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_eeg_ffilter_1(testCase)

spm('defaults','eeg');

fname = fullfile(spm('Dir'),'tests','data','OPM','test_opm.mat');
D     = spm_eeg_load(fname);
D = chantype(D,1:110,'MEG');
D.save();

% create 40Hz sinusoid
test = repmat(sin(2*pi*40*(0:.001:.999)),110,1);
D(:,:,1)= test;
D.save();

% try and remove sinusoid 
S=[];
S.D=D;
S.band='low';
S.freq=4;
fD = spm_eeg_ffilter(S);

% check that amplitude reduced by at least 100dB
Y = mean(std(D(:,900:1000,1),[],2));
res = mean(std(fD(:,900:1000,1),[],2));

dB = 20*log10(mean(Y./res));

testCase.verifyTrue(dB>100);
