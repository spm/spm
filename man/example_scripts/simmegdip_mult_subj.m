%% simulate data from multiple subjects. At the moment use the same dataset
%% template, 
%% can vary source std in amplitude and/or location across subjects

cd('D:\SPM MEEG 2010\Preprocessing demo\MEG');

Nsubj=10;
%% ORIGINAL FACES DATA FILE
origspmfilename='cdbespm8_SPM_CTF_MEG_example_faces1_3D.mat'
outpath='D:\Data\Scaling\simulated';
%% FOR SIMULATED NOISE- recording BW of 80Hz, white noise level 10ft/rtHz
noiselevel=50*1e-15;
BW=80; 
cross_subj_disp=5; %% sd of positions across subjects mm
cross_subj_sdamp_pc=100; %% %sd of amplitude variation across subjects

%% set frequency and amplitude of the two dipoles defined above 
%% for each condition (faces/scrambled) separately
    %cond1 cond2
dipfreq=[40 40;... %% dip 1
         20 20
         30 30 ];     %%% dip 2
      %cond1 cond2
dipamp=[1 0;...         % dip 1
        0.5 0
        1 0].*1e-1; % dip2 

    Ndips=size(dipfreq,1);
 substr=sprintf('simdata_%ddip_%dmm_%dpcamp',Ndips,cross_subj_disp,cross_subj_sdamp_pc);
%% define period over which dipoles are active
startf1=0.1; % (sec) start time
duration=0.3;% duration 




for s=1:Nsubj,
    dipamp_subj=dipamp+dipamp.*randn(size(dipamp))*cross_subj_sdamp_pc./100;

%% OUTPUT SIMULATED DATA FILE
spmfilename=[outpath filesep sprintf('%s%d',substr,s)];



%% FOR SIMULATED DATA
start_dipolepositions=[ 52, -29, 13; -52, -29, 13;-48 -22 48]; % in mni space
dipolepositions=start_dipolepositions+randn(size(start_dipolepositions)).*cross_subj_disp;
subj_data(s).dipole_positions=dipolepositions;


 
Ndips=size(dipolepositions,1);

%%% LOAD IN ORGINAL FACE DATA
D = spm_eeg_load(origspmfilename);
modality='MEG';
channel_labels = D.chanlabels(D.meegchannels(modality));

%%% GET LOCATION OF MESH Vertices in MNI and MEG/CTF space
allmeshvert_mni=D.inv{1}.mesh.tess_mni.vert;
allmeshvert_ctf=D.inv{1}.forward.mesh.vert;
allmeshfaces=D.inv{1}.forward.mesh.face;
allmeshnorms_ctf=spm_mesh_normals(struct('faces',allmeshfaces,'vertices',allmeshvert_ctf),true);
allmeshnorms_mni=spm_mesh_normals(struct('faces',allmeshfaces,'vertices',allmeshvert_mni),true);

%%% FORCE  A SINGE SPHERE HEAD MODEL- simpler to compare inversions
headmodels = {'Single Sphere', 'MEG Local Spheres', 'Single Shell'};
 D.inv{1}.forward(1).voltype=headmodels{1}; 
 D = spm_eeg_inv_forward(D); 
 grad = D.sensors('meg');
 vol = D.inv{1}.forward.vol;
 
 Ntrials=D.ntrials;

[condnames,ncond,trialtypes]=unique(D.conditions);
if length(condnames)~=size(dipamp_subj,2),
    error('number trial types should equal number of columns in dipamp');
end; % if


allavsignal=zeros(Ndips,length(D.time));
    
for dipind=1:Ndips,
cfg      = [];
cfg.vol  = vol;             
cfg.grad = grad;            

%% SIMULATE DIPOLES ON THE CORTICAL SURFACE


%% find nearest dipole to location specified
 [d meshind] = min(sum([allmeshvert_mni(:,1) - dipolepositions(dipind,1), ...
                             allmeshvert_mni(:,2) - dipolepositions(dipind,2), ...
                             allmeshvert_mni(:,3) - dipolepositions(dipind,3)].^2,2));
 meshdippos(dipind,:)=allmeshvert_ctf(meshind,:);
  cfg.dip.pos = meshdippos(dipind,:);
  t1=allmeshnorms_ctf(meshind,:); %% get dip orientation from mesh
  meshsourceind(dipind)=meshind;  %





cfg.dip.mom =t1';     % note, it should be transposed
cfg.ntrials =Ntrials;

endf1=duration+startf1; 

f1ind=intersect(find(D.time>startf1),find(D.time<=endf1));


for i=1:cfg.ntrials
    f1=dipfreq(dipind,trialtypes(i)); %% frequency depends on stim condition
    amp1=dipamp_subj(dipind,trialtypes(i));
    phase1=pi/2; 
    signal=zeros(1,length(D.time));
    signal(f1ind)=signal(f1ind)+amp1*sin((D.time(f1ind)-D.time(min(f1ind)))*f1*2*pi+phase1);
    cfg.dip.signal{i}=signal;
    allavsignal(dipind,:)=allavsignal(dipind,:)+signal;
end; % for i

cfg.triallength =max(D.time)-min(D.time);        % seconds
cfg.fsample = D.fsample;          % Hz
onesampletime=1/D.fsample;
cfg.absnoise=0;
cfg.relnoise=0;
cfg.channel=channel_labels;
raw1 = ft_dipolesimulation(cfg);
allraw(dipind)=raw1;
end; % for dipind

%% MERGE N DIPOLES INTO ONE DATASET
for dipind=1:Ndips,
    if dipind==1;
        dipsum=allraw(dipind).trial;
    else
        for i=1:Ntrials,
            newdata=cell2mat(allraw(dipind).trial(i));
            dataprev=cell2mat(dipsum(i));
            dataraw=newdata+dataprev;
            dipsum(i)=mat2cell(dataraw,size(newdata,1),size(newdata,2));
        end;
        
    end;
end;
allraw=raw1;
allraw.trial=dipsum;
clear dipsum;

    %% NOW ADD WHITE NOISE
    for i=1:Ntrials,
        dataraw=cell2mat(allraw.trial(i));
        channoise=randn(size(dataraw)).*noiselevel*sqrt(BW);
        dataraw=dataraw+channoise;
        allraw.trial(i)=mat2cell(dataraw, size(dataraw,1),size(dataraw,2));
        allraw.time(i)=mat2cell(D.time,1);
        epochdata=cell2mat(allraw.trial(i))';
    end;


avg1 = ft_timelockanalysis([], allraw);
plot(avg1.time, avg1.avg);  % plot the average timecourse

%%% now write out a new data set
D2=D;
D2=spm_eeg_ft2spm(allraw,spmfilename);
D2=sensors(D2,'MEG',raw1.grad);
D2=fiducials(D2,D.fiducials);
D2.inv=D.inv;
D2 = conditions(D2, [], D.conditions);
D2=coor2d(D2,'MEG',coor2d(D)); %% save projected 2d channel locations
D2.save;

figure;
h=plot(D2.time,allavsignal);
set(h(1),'Linestyle',':');
set(h,'LineWidth',4);
set(gca,'FontSize',18);
set(gcf,'color','w');

% %% write an averaged data set also
%   S=[]
%   S.D=D2;
%   S.robust=0;
%   D2av = spm_eeg_average(S);

end; % for s

