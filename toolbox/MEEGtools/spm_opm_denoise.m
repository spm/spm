function Dout = spm_opm_denoise(D,refD,derivative,gs,update,prefix)
% Denoise OPM data with specified regressors, optionally with derivatives
% and global signal.
% FORMAT Dout - spm_opm_denoise(D,refD,derivative,gs,update,prefix)
%
% D            - MEEG object
% refD         - array or MEEG object containing regresssors for denoising.
% derivative   - boolean indicating whether to use derivatives or not.
% gs           - boolean indicating whether to use Global Signal or not.
% update       - boolean indicating whether to create MEEG object.
% prefix       - string to prefix file path with if update is TRUE.
%__________________________________________________________________________
% Copyright (C) 2017-2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_opm_denoise.m 7414 2018-09-07 11:00:29Z spm $


% check to make sure dimensions match between reference and MEEG object
ndim = size(D);
ldim = size(refD);
if (~all(ndim(2:3)==ldim(2:3)))
    error('Dimensions of D and refD should be equal.');
end


% ref, sensors and residual objects 
ref = refD;
sensors = D;


% need to loop over Ntrials of size winSize
Ntrials = size(sensors,3);
winSize = size(sensors,2);
% only want to denoise MEG channels (not trigs or references)
megind=D.indchantype('MEG');
megres = zeros(size(sensors(megind,:,:)));

% loop (inefficient due to continued memory reallocation but ...)
for i =1:Ntrials
     % add a mean column to the reference regressors
     intercept = ones(winSize,1);
     reference = ref(:,:,i)';
     reference=[reference intercept];
    
     % optionally add derivatives
     if(derivative)
        drefer=diff(reference);
        drefer=[drefer(1,:); drefer];
        reference=[drefer reference];
     end
     
    % optionally add global signal 
    if(gs)
        trial = D(megind,:,i)';
        gsrefer =  mean(trial,2);
        reference=[gsrefer reference];
    end
    % reference is column major so transpose sensors
    beta = pinv(reference)*sensors(megind,:,i)';
    megres(:,:,i) = (sensors(megind,:,i)'- reference*beta)';
end

res=sensors(:,:,:); %% new file will only have MEG sensors changed
res(megind,:,:)=megres;

% make sure output has the singleton dimension if matrix supplied
if ((length(size(res)))==2)
    outsize = [size(res) 1];
else
    outsize = size(res);
end


if(update)
    % name, clone and fill with residual data 
    inFile = fnamedat(D);
    [a, b]=fileparts(inFile);
    outfile = fullfile(a,[prefix b '.dat']);
    Dout = clone(D,outfile,outsize);
    Dout(:,:,:) = res;
else
    Dout = res;
end
