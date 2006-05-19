function D = spm_eeg_artefact_TF(S)
% simple artefact detection
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%Stefan Kiebel, Rik Henson & James Kilner
% $Id: spm_eeg_artefact.m 265 2005-10-19 17:24:54Z james $


[Finter,Fgraph,CmdLine] = spm('FnUIsetup', 'EEG artefact setup',0);

try
    D = S.D;
catch
    D = spm_select(1, '.*\.mat$', 'Select EEG mat file');
end

P = spm_str_manip(D, 'H');

try
    D = spm_eeg_ldata(D);
catch
    error(sprintf('Trouble reading file %s', D));
end
try
    D2 = S.D2;
catch
    D2 = spm_select(1, '.*\.mat$', 'Select EEG mat file');
end

P = spm_str_manip(D, 'H');

try
    D2 = spm_eeg_ldata(D2);
catch
    error(sprintf('Trouble reading file %s', D));
end

Smoothing=round(Smoothing/1000*D.Radc);
MustDoWork = 1; % flag to indicate whether user already specified full artefact list



spm('Clear',Finter, Fgraph);

[Finter,Fgraph,CmdLine] = spm('FnUIsetup', 'EEG artefact setup',0);


D.events.reject = zeros(1, D.Nevents);

% cell vectors of channel-wise indices for thresholded events
D.channels.thresholded = cell(1, D.Nchannels);
index = [];
if MustDoWork 
    
    tloops = [1:D.Nchannels];
    tempsdata=zeros(D.Nchannels,D.Nfrequencies,D.Nsamples,length(D.events.Ntypes));
    
    for i = 1:D.events.Ntypes
        nbars=D.events.Ntypes*length(tloops);
        spm_progress_bar('Init', nbars, '2nd pass - robust averaging'); drawnow;
        if nbars > 100, Ibar = floor(linspace(1, nbars,100));
        else, Ibar = [1:nbars]; end
        
        trials = find(D.events.code == D.events.types(i));
        for nf=1:D.Nfrequencies
           for j = tloops %loop across electrodes	
         
                if ismember((i-1)*length(tloops)+j, Ibar)
                    spm_progress_bar('Set', (i-1)*length(tloops)+j);
                    drawnow;
                end
           
                tdata = squeeze(D.data(j, nf,:, trials));
                inds= find(tdata==0);  
                nt=0;
                for n=1:length(inds);
                    if inds(n)==1
                        nt=nt+1;
                        while inds(n+nt)==nt+n
                            nt=nt+1;
                        end  
                         tdata(inds(n))=tdata(inds(n)+nt);
                        
                    else
                        
                        tdata(inds(n))=tdata(inds(n)-1);
                    end
                end
                [B, bc] = spm_eeg_robust_averaget(log10(tdata),Weightingfunction,Smoothing);
                tempsdata(j,nf,:,i)=reshape(B,1,1,length(B),1);
                
             
            end
        end
        
        
    end
    
    spm_progress_bar('Clear');
    
    D.weights = allWf;
    
end % MustDoWork

% User-specified lists override any artefact classification
if D.thresholds.External_list
    D.events.reject(D.thresholds.out_list) = 1;
    D.events.reject(D.thresholds.in_list)  = 0;
end

% Save the data
copyfile(fullfile(D.path, D.fnamedat), fullfile(D.path, ['a' D.fnamedat]));
D.fnamedat = ['a' D.fnamedat];

spm_progress_bar('Clear');

D.data = [];
D.fname = ['a' D.fname];

if spm_matlab_version_chk('7') >= 0
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

spm('Pointer', 'Arrow');
