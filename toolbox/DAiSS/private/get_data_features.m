
function [trialfeatures,vedata]=get_data_features(flatpepdata,Nbins,Ntrials,weights,Tfull,datatype,featureind,regressout)

%%% put data in form for mv test
%% options are datatye
%% {'peakabs','peakreal','peakimag','peaksin','peakcos','sumpower'};
%% datatype=pwr power (re^2+comp^2)
%% datatpe=pwrsincos: power+sin and cos terms
%% datatype=4:
%% if weights==-1 just set up static variables depnding on datatype
%% 
flagnames={'peakabs','peakreal','peakimag','peaksin','peakcos','sumpower'};

if nargin==0, %% for menu
    trialfeatures=flagnames;
    return;
end;

if nargin<8
    regressout=[];
end;

persistent colflags


Nfeaturebands=numel(featureind); %% number of bands to take features from
%Ntrials=size(allepdata,1);
Nchans=size(flatpepdata,2);
%Nbins=size(allepdata,3);


if weights==-1,
    colflags=zeros(1,numel(flagnames));
    disp('using the following features per band');
    for f=1:numel(flagnames),
        colflags(f)=~isempty(findstr(datatype,flagnames{f}));
        if colflags(f),
            disp(sprintf('%s',flagnames{f}));
        end;
    end;

trialfeatures=[];
return;
end; % if weights=-1

fperbandind=find(colflags);
Nfperband=length(fperbandind);

trialfeatures=zeros(Ntrials,Nfeaturebands*Nfperband);


vedata=flatpepdata*weights';
vedata=vedata-mean(vedata);
if ~isempty(regressout),
    %keyboard;
     [B,BINT,residuals,rint,stats] = regress(vedata,[regressout ones(size(regressout)).*std(regressout)]);
     vedata=residuals;
     if stats(1)==1, %% r squared is equal to one, i.e. 100% variance explained
         disp('replacing perfectly predicted data with random sequence')
         vedata=randn(size(vedata));
     end;
end;
trialvedata=(reshape(vedata,Nbins,Ntrials))';

trialvedata=trialvedata-repmat(mean(trialvedata,2),1,Nbins); %% remove mean
trialvedata=trialvedata.*repmat(spm_hanning(Nbins)',Ntrials,1); %% window

features=zeros(Ntrials,Nfperband);
for f=1:Nfeaturebands,
    Tband=Tfull(:,featureind{f}); % filter to this band
    
    dcttrialdata=trialvedata*Tband; %% frequency representation
    ftrialdata=dcttrialdata*Tband'; %% filter data
    
    Ff=fft(ftrialdata');
    meanpower=mean(abs(Ff')); %% get bin with highest power and use the phase from this
    [maxval,maxbin]=max(meanpower(round(1:length(meanpower)/2)+1)); %% only take 1st half
    
    features(1:Ntrials,1)=abs(Ff(maxbin,:))'; %% could use mean or peak
    features(1:Ntrials,2)=real(Ff(maxbin,:))';
    features(1:Ntrials,3)=imag(Ff(maxbin,:))';
    features(1:Ntrials,4)=sin(angle(Ff(maxbin,:)))';
    features(1:Ntrials,5)=cos(angle(Ff(maxbin,:)))';
    
    features(1:Ntrials,6)=sum(Ff.*conj(Ff))';
    trialfeatures(1:Ntrials,(f-1)*Nfperband+1:f*Nfperband)=features(1:Ntrials,fperbandind);
    
end; % for f






