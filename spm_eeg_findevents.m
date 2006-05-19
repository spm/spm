function D = spm_eeg_findevents(S)
% function used for epoching continuous EEG/MEG data
% FORMAT D = spm_eeg_epochs(S)
% 
% S		    - optional input struct
% (optional) fields of S:
% D			- filename of EEG mat-file with continuous data
% events    - struct with various entries:
%    start     - pre-stimulus start of epoch [ms]
%    stop	   - post-stimulus end of epoch [ms]
%    types	   - events to extract (vector of event types)
%    Inewlist  - switch (0/1) to have new list of event codes
%    Ec        - list of new event codes
%
% Output:
% D			- EEG data struct (also written to files)
%_______________________________________________________________________
%
% spm_eeg_epochs extracts single trials from continuous EEG/MEG data. The
% length of an epoch is determined by the samples before and after stimulus
% presentation. One can limit the extracted trials to specific trial types.
% Also, it is possible to re-number trial types, see above. 
% Note that epoching includes a baseline correction of each single trial,
% i.e. a subtraction of the average pre-stimulus average from all time
% points.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_epochs.m 213 2005-08-22 12:43:29Z stefan $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG epoching setup',0);



try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
end

P = spm_str_manip(D, 'H');

try
	D = spm_eeg_ldata(D);
catch    
	error(sprintf('Trouble reading file %s', D));
end

if ~isfield(D, 'events')
    D.events = [];
end
	



% transform ms to samples
stops = ceil(2000*D.Radc/1000);

% Allocate space for epoched data


d = zeros(D.Nchannels, stops);


stims=floor(D.Nsamples/stops);


[Cel, Cind, x, y] = spm_eeg_locate_channels(D, 128, 1);
d1=[];

figure(2)


text(Cel(:,1),Cel(:,2),[D.channels.name(:)])
axis([-10 150 -10 150])

but =1;

channels=zeros(1,275);
while but == 1
    [x1,y1,but]=ginput(1);
    tmp=((Cel(:,1)-x1).^2+(Cel(:,2)-y1).^2).^0.5;
    [junk,elec]=sort(tmp);
    if but ==1
        if channels(elec(1))==0
            channels(elec(1))=1;
        else
            channels(elec(1))=0;
        end
    end
    clf
    text(Cel(find(channels==0),1),Cel(find(channels==0),2),[D.channels.name(find(channels==0))])
    axis([-10 150 -10 150])
end
inds=find(channels==1);
ev=10;
D.events.time= [];
D.events.code= [];
clf
for jk=1:stims
    
    
    d = D.data(inds, (stops*(jk-1))+1 :stops*jk, 1);
    
    % baseline subtraction
    
    d = d - repmat(mean(d(:,:),2), 1, stops);
    figure(2)
    plot(d')
    title(['trial ' num2str(jk) ' out of ' num2str(stims) ' trials'])
    
    [x,y,but]=ginput(1);
    
    while (but(1)~=3)
        
        D.events.time= [D.events.time, (stops*(jk-1))+x];
        D.events.code= [D.events.code, ev];
      
        hold on
        plot([x,x],[min(min(d)),max(max(d))],'k')
        ev=ev+1;
        [x,y,but]=ginput(1);
      
    end
    clf
end
        

            
 



if spm_matlab_version_chk('7') >= 0
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

