function D = spm_opm_epoch_trigger(S)
% Epoch M/EEG data based on triggers (channel type is TRIG) in the file.
% FORMAT D = spm_opm_epoch_trigger(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                    - Default: no Default
%   S.timewin       - n x 2 matrix where n is the        - Default: replicates the first two numbers for each condition
%                     numer of conditions and the 
%                     2 numbers are the time around
%                     the trigger in ms.
%   S.condLabels    - n x 1 cell containing condition    -Default: Cond N
%                     labels
% Output:
%  D           - epoched MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_opm_epoch_trigger.m 7419 2018-09-14 13:17:03Z tim $


%-Set Defaults
%--------------------------------------------------------------------------
if ~isfield(S, 'D'),       error('D is required'); end
cTypes= chantype(S.D);
trigInds= strmatch('TRIG',cTypes);
nTrigs =length(trigInds);
trigs=S.D(trigInds,:);

if ~isfield(S, 'condLabels')
    args =[];
    args.base='Cond';
    args.n=nTrigs;
    S.condLabels=spm_create_labels(args);
end

if(size(S.timewin,1)<nTrigs)
    S.timewin= repmat(S.timewin,nTrigs,1);
end

%-Create trl and cond matrices
%--------------------------------------------------------------------------
trl=zeros(0,3);
cond ={};
nevents=[];
for i=1:nTrigs
    
    tChan=trigs(i,:);                                  % Get ith trigger
    thresh=std(tChan)*2;                               % Threshold it
    
    evSamples=find(diff(tChan>thresh)==1)+1;           % Find 1st sample exceeding threshold
    offsetTime=S.timewin(i,1)/1000;                    % Offset in seconds
    offsetSamples=round(offsetTime.*S.D.fsample);      % Offset in samples
    durationTime=diff(S.timewin(i,:))/1000;            % Duration in seconds
    durationSamples=round(durationTime.*S.D.fsample);  % Duration in samples
    
    begSample=evSamples+offsetSamples;                 % Get first sample
    endSample=begSample+durationSamples;               % Get end sample
    offset=offsetSamples.*ones(size(begSample));       % Replicate offset accross events
    
    trlTemp=round([begSample'  endSample'  offset']);  % Construct temporary trial matrix
    nevents(i)=size(trlTemp,1);                        % Compute number of events per Condition
    trl = [trl;trlTemp];                               % Combine trl matrices accross conditions
    
    condTemp =repmat({S.condLabels{i}},nevents(i),1);  % Replicate codition lable accross events
    cond = {cond{:,:},condTemp{:,1}}';                 % Combine condition labels accross conditions
    
end


% Actually do the epoching now
%--------------------------------------------------------------------------
args = [];
args.D = S.D;
args.trl = trl;
args.conditionlabels =cond;
args.bc = 0;
args.prefix = 'e_';
D = spm_eeg_epochs(args);
