function [mfD,Yinds] = spm_opm_hfc(S)
% remove interference that behaves as if it was from a harmonic (magnetic) field
% FORMAT D = spm_opm_hfc(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                                - Default: no Default
%   S.L             - Spherical harmonic order (1=homogenous field)  - Default: 1
%   S.usebadchans   - logical to correct channels marked as bad      - Default: 0
%   S.chunkSize     - max memory usage(for large datasets)           - Default 512(MB)
%   S.badChanThresh - threshold (std) to identify odd channels       - Default 50 (pT)
%   S.balance       - logical to update forward model                - Default 1
%   S.prefix        - prefix to filename                             - Default 'h'
% Output:
%   D               - denoised MEEG object (also written to disk)
%   Yinds           - the indices of filtered channels
%__________________________________________________________________________

% Tim Tierney
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging


%-Set default values
%--------------------------------------------------------------------------
errorMsg = 'an MEEG object must be supplied.';
if ~isfield(S, 'D'),             error(errorMsg); end
if ~isfield(S, 'usebadchans'),   S.usebadchans = 0; end
if ~isfield(S, 'chunkSize'),     S.chunkSize = 512; end
if ~isfield(S, 'badChanThresh'), S.badChanThresh = 50; end
if ~isfield(S, 'balance'),       S.balance = 1; end
if ~isfield(S, 'L'),             S.L = 1; end
if ~isfield(S, 'prefix'),        S.prefix = 'h'; end

%-Find usable channels
%--------------------------------------------------------------------------

s = sensors(S.D,'MEG');
if isempty(s)==1
    error('Could not find sensor positions')
end

chaninds = indchantype(S.D,'MEG');
if(S.usebadchans)
  usedLabs= chanlabels(S.D,chaninds);
else
  badinds = badchannels(S.D);
  usedinds = setdiff(chaninds,badinds);
  usedLabs= chanlabels(S.D,usedinds);
end

usedLabs = intersect(usedLabs,s.label);
[~,sinds] = spm_match_str(usedLabs,s.label);

%-Define harmonic Basis Set
%--------------------------------------------------------------------------
args=[];
args.o= s.coilori(sinds,:);
args.v = s.coilpos(sinds,:);
args.li = S.L;
X = spm_opm_vslm(args);
fprintf('%-40s: %30s\n','Created Design Matrix',spm('time'));

%-Compute projector
%--------------------------------------------------------------------------
M = eye(size(X,1))-X*pinv(X);

%-Get Data indices
%--------------------------------------------------------------------------
Yinds = indchannel(S.D,usedLabs);

if (size(Yinds,1)~=size(X,1))
    error('data size ~= number of sensors with orientation information');
else
    
%-create ouput dataset object
%--------------------------------------------------------------------------
fprintf('Creating output dataset\n'); 
outname = fullfile(path(S.D),[S.prefix fname(S.D)]);
mfD = clone(S.D,outname);
mfD.save();

%- Work out chunk size
%--------------------------------------------------------------------------
chunkSamples= round(S.chunkSize/(8*size(S.D,1))*1e6);
begs=1:chunkSamples:size(S.D,2);
ends = (begs+chunkSamples-1);
if(ends(end)>size(S.D,2))
    ends(end)= size(S.D,2);
end

%-Run on channels needing correction
%--------------------------------------------------------------------------
trvar= zeros(length(Yinds),size(S.D,3));
fprintf('%-40s: %30s\n','Processing Data',spm('time'));

for j=1:size(S.D,3)
    mk= S.D(Yinds,1,j);
    sk = zeros(length(Yinds),1);
    count=1;
    
    for i =1:length(begs)
        inds = begs(i):ends(i);
        
        % mfD(Yinds,inds,j)=M*S.D(Yinds,inds,j) is slow (disk read)
        out = S.D(:,inds,j);
        Y=out(Yinds,:);
        out(Yinds,:)=M*Y;
        mfD(:,inds,j)=out;
        
        % accurate running variance for identifying odd channels
        % (https://www.johndcook.com/blog/standard_deviation/)
        for l = 1:length(inds)
            xk = out(Yinds,l);
            mkprev = mk;
            mk = mkprev +(xk-mkprev)/count;
            sk=sk+(xk-mkprev).*(xk-mk) ;
            count=count+1;
        end
        
    end
    trvar(:,j)=sk/(count-1);
end

%-Update forward modelling information
%--------------------------------------------------------------------------
if (S.balance)
    fprintf('%-40s: %30s\n','Updating Sensor Information',spm('time'));
    grad = mfD.sensors('MEG');
    tmpTra= eye(size(grad.coilori,1));
    tmpTra(sinds,sinds)=M;
    grad.tra                = tmpTra*grad.tra;
    grad.balance.previous   = grad.balance.current;
    grad.balance.current    = 'hfc';
    mfD = sensors(mfD,'MEG',grad);
    % Check if any information in D.inv needs updating.
    % TODO: Update to support multiple invs/forwards/modalities
    if isfield(mfD,'inv')
        if isfield(mfD.inv{1},'gainmat')
            fprintf(['Clearing current forward model, please recalculate '...
                'with spm_eeg_lgainmat\n']);
            mfD.inv{1} = rmfield(mfD.inv{1},'gainmat');
        end
        if isfield(mfD.inv{1},'datareg')
            mfD.inv{1}.datareg.sensors = grad;
        end
        if isfield(mfD.inv{1},'forward')
            voltype = mfD.inv{1}.forward.voltype;
            mfD.inv{1}.forward = [];
            mfD.inv{1}.forward.voltype = voltype;
            mfD = spm_eeg_inv_forward(mfD,1);
        end
    end
    mfD.save();
end

%-Odd Channel Check
%--------------------------------------------------------------------------
fprintf('Checking for unusual channels\n');
SD = mean(sqrt(trvar),2)*1e-3;
for i = 1:length(SD)
    index= indchannel(S.D,usedLabs{i});
    if(SD(i)>S.badChanThresh)
        fprintf(['Residual on channel ' num2str(index) ', '...
            usedLabs{i} ': %3.2f pT\n'], SD(i));
    end
end

%-Complete
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Completed',spm('time'));

end
