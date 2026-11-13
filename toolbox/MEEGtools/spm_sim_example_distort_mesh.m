% Clean up
clear all;
close all;

% Set defaults SPM settings
spm('defaults','eeg');

% Set results folder
%resultsdir='D:\matlab\jimmymatfiles\'; %% user specified.
distortdir='D:\spm_distort'
addpath(distortdir)
ar=load('righthemind.mat');
al=load('lefthemind.mat');

resultsdir = fullfile(distortdir, 'sim_results');  % standard within current folder, but user should adapt


% Flexibly create results folder if missing
if not(isfolder(resultsdir))
    mkdir(resultsdir)
end

% Main parameter
REDOCALC=1; %% should =1 first time code is run or if want to simulate again (Or calc new meshes etc)

%% After that = 0 just loads in results.

if ~exist(resultsdir),
    error('need to set up a directory to hold the results !')
end;

fname = fullfile(spm('Dir'),'tests','data','meg_face','cdbespm8_SPM_CTF_MEG_example_faces1_3D.mat');

if ~exist(fname),
    error('Need to download SPM test data')
end;

D = spm_eeg_load(fname);

% Fix paths to make sure to use the canonical templates included in one's
% own SPM version
D = spm_setmeshpaths(D,fullfile(spm('Dir'),'canonical'))

invmethods={'EBB','IID'}; %% methods to be used

%% Simulate a source
prefix='sim';
freq=30; %% frequency of sinusoid
scale=[0.1]; % how much to scale simulated source by

HeadModel='Single Shell'; % forward model to use
peakdistthresh=20; %% min distance between local maxima to be deemed different

%% Now make some distorted brians.
Npoints=17; %% number of surfaces
Zrange=3; %% amount of distortion (as Z score)
Nseed=8; %% number of random seeds to use
Zrange=[-Zrange,Zrange];
trajvals=linspace(Zrange(1),Zrange(2),Npoints);
DistortIndices=[8 100]; %% only distort higher order components

Crossval=0; %% do cross validation or not (takes ~10* longer)

simstr={'1motor'}; %% simulation scenario


Jpk=[];
Fieldmap=[];
Lf=[];

if Crossval,
    Nblocks=10; pctest=10;
else
    Nblocks=1; pctest=0; %% not testing cross validation, but could do..
end;


% Loop over sims
if REDOCALC

    for simind=1:numel(simstr), %% was on simind=12
        simpos=simstr{simind}
        switch (simpos),
            case '1motor',
                Nsources=1;
                mnicoord=[-40, -18, 48]; %% motor cortex approx
                mult_simsignal=sin(freq*[1:D.nsamples]*2*pi/D.nsamples).*scale; %% single sinusoid cycle


                %% Can add other scenarios here
            otherwise
                error('case not found')

        end; % case simpos

        %% Find a peak in the simulated signal at around the central cycle (only used for visualization in figure)
        l1=length(mult_simsignal(1,:));
        l2=round(freq*2*pi/D.nsamples)*4;
        [dum,peakl2]=max(mult_simsignal(1,round(l1/2)-l2:round(l1/2)+l2))
        pkind=peakl2+round(l1/2);


        [Dnew,meshsourceind]=spm_eeg_simulate({D},prefix,mnicoord,mult_simsignal);

        PCAtemplate=fullfile(spm('Dir'),'tpm','shp','Template_0.nii');

        F_vals=[];
        R2=[];
        maxsumJ2=[];

        derr=zeros(Nseed,Npoints);

        for RandSeed=1:Nseed,

            matlabbatch=[];
            matlabbatch{1}.spm.meeg.source.eeg_shp_distort.D = {Dnew.fullfile};
            matlabbatch{1}.spm.meeg.source.eeg_shp_distort.val = 1;
            matlabbatch{1}.spm.meeg.source.eeg_shp_distort.PCAtemplate = {PCAtemplate};
            matlabbatch{1}.spm.meeg.source.eeg_shp_distort.TemplateRedo = 'No';
            matlabbatch{1}.spm.meeg.source.eeg_shp_distort.DistortIndices = DistortIndices;
            matlabbatch{1}.spm.meeg.source.eeg_shp_distort.Npoints = Npoints; %% only 4 points on trajectory
            matlabbatch{1}.spm.meeg.source.eeg_shp_distort.Zrange = Zrange; %% move from minimal to maximal distortion
            matlabbatch{1}.spm.meeg.source.eeg_shp_distort.Zlimit = 6.5;
            matlabbatch{1}.spm.meeg.source.eeg_shp_distort.Zcentre = 'sub'; %% keep it subject specific
            matlabbatch{1}.spm.meeg.source.eeg_shp_distort.RandSeed = RandSeed; %% random seed defines trajectory

            [aM,b]=spm_jobman('run', matlabbatch);

            %% Now compute distances of brians to brain.
            M0=gifti(aM{1}.brain); %% get original cortex
            for i=1:numel(aM{1}.brians),
                M1=gifti(aM{1}.brians(i));
                derr(RandSeed,i)=mean(sqrt(dot((M1.vertices-M0.vertices)',(M1.vertices-M0.vertices)')));
            end;


            %% Now make lead fields for each surface (looks in directory seedXXX as specified by RandomSeed above)

            LFsubdir=sprintf('seed%03d',RandSeed);
            matlabbatch=[];
            matlabbatch{1}.spm.meeg.source.eeg_shp_gainmat.D = {Dnew.fullfile};
            matlabbatch{1}.spm.meeg.source.eeg_shp_gainmat.val = 1;
            matlabbatch{1}.spm.meeg.source.eeg_shp_gainmat.LFheadmodel = HeadModel;
            matlabbatch{1}.spm.meeg.source.eeg_shp_gainmat.LFRedo = 'Yes';
            matlabbatch{1}.spm.meeg.source.eeg_shp_gainmat.LFsubdir = LFsubdir;

            [aL,b]=spm_jobman('run', matlabbatch);


            %% Now quantify how much the lead fields for the distorted surfaces differ from the original
            G=load(aL{1}.gainmat{1}); %% original lead fields
            L0=G.G;
            for i=1:Npoints+1,
                G=load(aL{1}.gainmat{i});
                L1=G.G;
                diffL(i)=mean(dot(L1-L0,L1-L0));
            end; % for i


            %% Now run an inversion for the new brians and the original brain
            val=1;
            fprintf('\n Using these gainmats to make generic U matrix\n')
            gainmatfiles=[];
            for f3=1:numel(aL{1}.gainmat)
                gainmatfiles=strvcat(gainmatfiles,aL{1}.gainmat{f3});
            end;


            % Setup spatial modes for cross validation and to ensure same modes used
            % across the group of inversions (not biased to first or last etc)
            idstr='testsim';

            [a1 b1 c1]=spm_fileparts(gainmatfiles(2,:));
            spatialmodesname=fullfile(a1, sprintf('Grpmod_%s_seed%03dNb%d.mat',idstr,RandSeed,Npoints));
            [spatialmodesname,Nmodes,pctest]=spm_eeg_inv_prep_modes_xval(Dnew.fullfile, [], spatialmodesname, Nblocks, pctest,gainmatfiles);


            %% Now invert these data using this headmodel

            Ngainmat=size(gainmatfiles,1);

            foi=[0 freq*2]; %% frequency range
            woi=[Dnew.time(1) Dnew.time(end)].*1000; %% time window in ms
            patch_size=0; % should correspond to simulated dipole extent above (defaults to 0)
            n_temp_modes=[]; %% number of temporal modes to use

            for f=1:Ngainmat,
                starttime=tic;
                D=spm_eeg_load(Dnew.fullfile);

                %% Load in a dataset but swap in different gainmat (lead field) files
                gainmatfile=deblank(gainmatfiles(f,:));
                [gpath,gname,gext]=spm_fileparts(gainmatfile);
                copyfile(gainmatfile,[D.path filesep gname gext]); % have to put gainmat in spmfile directory
                D.inv{val}.gainmat=[gname gext]; %%  just change name of gainmat in file
                D.save;

                %% Now run the inversion
                matlabbatch={};
                batch_idx=1;
                for i1=1:numel(invmethods),

                    % Source reconstruction
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.D = {Dnew.fullfile};
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.val = 1;
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.whatconditions.all = 1;
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.invfunc = 'Classic';
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.invtype = invmethods{i1}; %;
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.woi = woi;
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.foi = foi;
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.hanning = 0;
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.randpatch.npatches = 512;
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.randpatch.niter = 1;
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.patchfwhm = patch_size; %% NB A fiddle here- need to properly quantify
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.mselect = 0;
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.nsmodes = Nmodes;
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.umodes = {spatialmodesname};
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.ntmodes = n_temp_modes;
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.priors.priorsmask = {''};
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.priors.space = 1;
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.restrict.locs = zeros(0, 3);
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.restrict.radius = 32;
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.outinv = '';
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.modality = {'All'};
                    matlabbatch{batch_idx}.spm.meeg.source.invertiter.crossval = [pctest Nblocks];
                    spm_jobman('run', matlabbatch);

                    % Load in reconstruction. Get metrics of inversion
                    Drecon=spm_eeg_load(Dnew.fullfile);
                    F_vals(f,i1,RandSeed)=Drecon.inv{1}.inverse.F;

                    R2(f,i1,RandSeed)=Drecon.inv{1}.inverse.R2;
                    VE(f,i1,RandSeed)=Drecon.inv{1}.inverse.VE;
                    crosserr(f,i1,RandSeed)=mean(Drecon.inv{1}.inverse.crosserr(:));
                    T=Drecon.inv{1}.inverse.T; % temporal modes
                    allJ=Drecon.inv{1}.inverse.J{1}*T';
                    if RandSeed==1, %% data for figure;
                        L=Drecon.inv{1}.inverse.L;
                        Jpk(f,i1,:)=allJ(:,pkind); %% current distribution at peaktime
                        Fieldmap(f,i1,:)=L*allJ(:,pkind); %% predicted field due to the inversion
                        Lf(f,:,:)=L;
                        c2d=Drecon.coor2D;
                        chanlabel=Drecon.chanlabels;
                    end;
                    %% Now get classical distances between simulated and actual sources

                    meshpwr=sum(allJ'.^2); %% total power across mesh
                    [maxsumJ2(f,i1,RandSeed),ind]=max(meshpwr); %%
                    Mmni=Drecon.inv{1}.mesh.tess_mni;
                    M=[];
                    M.faces=Mmni.face;
                    M.vertices=Mmni.vert;
                    L = spm_mesh_get_lm(M,meshpwr'); %% get local power maxima over mesh
                    [maxvals,localmaxinds]=sort(meshpwr(L),'descend'); %% sort in descending order
                    Lsort=L(localmaxinds);


                    useind=Lsort(1); %% for 1st source, this is the global maximum index
                    allmni=M.vertices(useind,:); %% for 1st source, this is the global maximum location in mni space
                    allvals=maxvals(1);


                    count=1;

                    for f2=1:Nsources-1, % for more than 1 source
                        dx=repmat(M.vertices(Lsort(f2+1),:),size(allmni,1),1)-allmni;
                        alldist=sqrt(dot(dx',dx')); %% get distance from this local max to last global maximum

                        if min(alldist)>peakdistthresh, %% if it is far enough from last maximum noted
                            count=count+1;
                            useind(count)=Lsort(f2+1); %% record index
                            allmni(count,:)=M.vertices(useind(count),:); %% record position
                            allvals(count,:)=maxvals(f2+1); %% and peak value
                        end; % if

                    end; % for f2
                    d0=[]; %% now look how close the simulated locations (allmni) were to the local maxima
                    for s=1:Nsources,
                        d1=(repmat(mnicoord(s,:),size(allmni,1),1)-allmni)';
                        d0(s)=min(sqrt(dot(d1,d1)));
                    end;

                    %% Record some distances between simulated and reconstructed
                    Mindist(f,i1,RandSeed)=min(d0);
                    Maxdist(f,i1,RandSeed)=max(d0);
                    Meandist(f,i1,RandSeed)=mean(d0);

                end; % for invmethods


                %% Clean up
                fprintf('\n Done %d of %d',f,Ngainmat);
                fprintf('\n Deleting Gainmat')
                delete([D.path filesep gname gext]); %% remove copied gainmat file to save disk space
                D.inv{val}.gainmat=''; %%  just change name of gainmat in spmfile
                D.save;

            end; % for f=1:Ngainmat (over surrogates)



            save([resultsdir filesep sprintf('summary_%d_cond_CTF%s_seed%d_scale%03d.mat',Crossval,simpos,RandSeed,round(scale*1000))])

            close all;
        end; % for RandSeed
    end % for simind
else %if REDOCALC=0
    RandSeed=Nseed;
    simpos=simstr{end};
    load([resultsdir filesep sprintf('summary_%d_cond_CTF%s_seed%d_scale%03d.mat',Crossval,simpos,RandSeed,round(scale*1000))])
end; % if REDOCALC

%% Plot preferences

mstr='do'; % markers for inversions
colstr='br'; %colours for inversions
col_dist_str='mc'; %% for distorted surfaces
lspec_dist='-:'; %% line spec for true vs. distorted
plotinv=[1,2]; %% methods to plot

disterr=[0 mean(derr)]; %average distance error, including true brain at start
simind=1;

simpos=simstr{simind};

usesurfind=[1,Npoints]; %% zero, max distortion.

for invind=1:numel(invmethods),
    plotinvmeth=invmethods(invind);



    c2d=Drecon.coor2D;

    chanlabel=Drecon.chanlabels;

    megind=Drecon.indchantype('MEG')
    Y=Drecon.inv{1}.inverse.Y;
    T=Drecon.inv{1}.inverse.T;
    YT=Y*T';

    meandata=YT(:,pkind)
    empmean=mean(Drecon(megind,pkind,:),3);


    M0=gifti(aM{1}.brain);
    d1=zeros(Ngainmat,1);

    M0i=spm_mesh_inflate(M0);
    allplotJ=[];
    for fs=1:length(usesurfind),
        f=usesurfind(fs);
        if f>1,
            Mmesh{f}=gifti(aM{1}.brians{f-1})
        else
            Mmesh{f}=M0;
        end;
        d0=Mmesh{f}.vertices-M0.vertices;
        d1(f)=mean(sqrt(dot(d0',d0')));

        XYZ=Drecon.inv{1}.mesh.tess_mni.vert;

        J=squeeze(Jpk(f,invind,:));

      
         plotJ=abs(J);
         plotJ(find(plotJ<max(plotJ/20)))=0;
         allplotJ(fs,:)=plotJ;

       [ZI]=spm_eeg_plotScalpData(squeeze(Fieldmap(fs,invind,:)),c2d,chanlabel);
        allPredmap(fs,:,:)=ZI;
        title('Predicted field')
        Lstartpoint=squeeze(Lf(f,:,meshsourceind));
        [ZI]=spm_eeg_plotScalpData(Lstartpoint',c2d,chanlabel);
        title('Fieldmap of original vertex')
        allFmap(fs,:)=squeeze(Fieldmap(f,invind,:));
        allZFmap(fs,:,:)=ZI;

    end; % for fs

    figure ; hold on;
    
    maxmap=max(max(abs(allPredmap(1,:,:))))
    mstep=2*maxmap/10;
    clines=[-maxmap+mstep:mstep:maxmap-mstep]
    for fs=1:length(usesurfind),
        [c,h]=contour(squeeze(allPredmap(fs,:,:)),clines,lspec_dist(fs)); %% J*L
        set(h,'Linewidth',4)

        if fs==2,
            set(h,'Linewidth',2)
            set(h,'Linecolor',col_dist_str(fs))
        end;
        

    end;
    fprintf('\nPc variance between distorted and undistorted sensor estimates for %s',invmethods{invind})
    corrcoef(Fieldmap(usesurfind(1),invind,:),Fieldmap(usesurfind(2),invind,:))
    fprintf('\nPc variance between undistorted sensor estimates for %s and %s',invmethods{1},invmethods{2})
    corrcoef(Fieldmap(usesurfind(1),1,:),Fieldmap(usesurfind(1),2,:))
    fprintf('\nPc variance between distorted sensor estimates for %s and %s',invmethods{1},invmethods{2})
    corrcoef(Fieldmap(usesurfind(2),1,:),Fieldmap(usesurfind(2),2,:))

    set(gca,'PlotBoxAspectRatio',[0.9 1 1])
    axis([10  80 10 100])
    h=colorbar;
    ylabel(h,'Predicted field map (fT)')
    set(gca,'Fontsize',18)
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    title([invmethods{invind}])
    figure ; hold on;
    
    maxmap=max(max(abs(allZFmap(1,:,:))))
    mstep=2*maxmap/10;
    clines=[-maxmap+mstep:mstep:maxmap-mstep]
    for fs=1:length(usesurfind),
        [c,h]=contour(squeeze(allZFmap(fs,:,:)),clines,lspec_dist(fs))
        set(h,'Linewidth',4)
        if fs==2,
            set(h,'Linewidth',2)
            set(h,'Linecolor',colstr(fs))
        end;

    end;
    h=colorbar;
    axis([10  80 10 100])

    set(gca,'PlotBoxAspectRatio',[0.9 1 1])
    ylabel(h,'Field (fT) per nAm')
    set(gca,'Fontsize',18)
    set(gca,'XTick',[])
    set(gca,'YTick',[])


    figure; hold on;
    slicedim=2;
    simmeshpos=M0.vertices(meshsourceind,:);
    sliceind=simmeshpos(:,slicedim);
    viewind=zeros(1,3);
    viewind(slicedim)=1;
    sthick=3;
    approxcentre=mean(M0.vertices);
    tancomp=[];
    slview=setdiff([1:3],slicedim);
    useind=intersect(find(M0.vertices(:,slicedim)<sliceind+sthick),find(M0.vertices(:,slicedim)>sliceind-sthick))
    C=zeros(size(M0.vertices,1),1);
    C(useind)=1;

    
    %% PLOTS DISTORTED ANATOMY
    leftsliceind=al.orderedsliceind;
    rightsliceind=ar.orderedsliceind;
    allind=[leftsliceind ;rightsliceind];

    for fs=1:length(usesurfind),
        f=usesurfind(fs);
        [Norm_vert] = spm_mesh_normals(Mmesh{f});
        simpos=Mmesh{f}.vertices(meshsourceind,:);
        simori=Norm_vert(meshsourceind,:);
        %pos=al.pos;
        pos=Mmesh{f}.vertices;
        h=plot(pos(allind,1),pos(allind,3),col_dist_str(fs));
        set(h,'Linewidth',3)
 
        hold on;
        qscale=6;
        if fs==1,
        h1=quiver(simpos(:,1),simpos(:,3),...
            simori(:,1),simori(:,3),qscale);
        set(h1,'Color',colstr(invind))
        set(h1,'Linewidth',4)
        end;

    end; % for fs
    set(gca,'xdir','reverse'); %% consisent with paper

    %% PLOTS DISTORTED ANATOMY WITH CURRENT MAG

    for fs=1:length(usesurfind),
        figure; hold on;
        f=usesurfind(fs);
        nv=spm_mesh_normals(Mmesh{f});
        nvslice=nv(find(C),:);
        MS = spm_mesh_split(Mmesh{f}, C);

        plotJslice=abs(allplotJ(fs,find(C)));
        plotJslice=plotJslice./max(plotJslice)
        highind=find(plotJslice>0.5);
        h=plot(pos(allind,1),pos(allind,3),col_dist_str(fs));
        set(h,'Linewidth',3)

        
        for g=1:length(highind),
            Jscale=10*plotJslice(highind(g));
            arrowscale=0.6;
            hq=quiver(MS.vertices(highind(g),1),MS.vertices(highind(g),3),...
                nvslice(highind(g),1)*Jscale,nvslice(highind(g),3)*Jscale,arrowscale);
            set(hq,'Linewidth',3);
            set(hq,'color',colstr(invind))
            
            axis([-70 -2 10 78])
            set(gca,'xdir','reverse')
            set(gca,'visible','off')
        end;
    end; % for f
end; % for invind

tmpFvals=F_vals(:,plotinv,1)
figure;
sFvals=[];
for f=1:length(plotinv),

    [sFvals(:,f),sdist]=simplify_dist_metric(tmpFvals(:,f)',disterr)
end;
hold on;

for f=1:length(plotinv)

    h=plot(sdist,sFvals(:,f),colstr(f))
    [val(f),ind(f)]=max(sFvals(:,f));
    set(h,'Linewidth',4)
    hold on;

end;

legstr=invmethods(plotinv);
for f=1:length(plotinv)
    h2=plot(sdist(ind(f)),val(f),[colstr(f) mstr(f)])
    set(h2,'MarkerSize',15)
    legstr{end+1}=['Peak ' invmethods{plotinv(f)}];
end;
set(gca,'FontSize',18)
ylabel('Free Energy')
xlabel('Distortion (mm)')
legend(legstr)




medists=[];
sedists=[]


figure;

hold on;
alledists=[];


for i1d=1:length(plotinv),
    i1=plotinv(i1d);
    [tF_vals,sdist]=simplify_dist_metric(squeeze(F_vals(:,i1,:))',disterr);
    h=plot(sdist,squeeze(tF_vals),colstr(i1d));
    set(h,'Linewidth',1)
    [maxval,maxinds]=max(squeeze(tF_vals)');

    h01=plot(sdist(maxinds(1)),maxval(1),[colstr(i1d) mstr(i1d)])
    set(h01,'Markersize',20)


    h0=plot(sdist(maxinds),maxval,[colstr(i1d) mstr(i1d)])
    set(h0,'Markersize',20)

    axis([0 max(sdist) -Inf Inf]);
    xlabel('Distortion (mm)')
    set(h(1),'Linewidth',2)
    ylabel('Free Energy')

    edists=sdist(maxinds);
    alledists(i1d,simind,:)=edists;
    medists(i1d,simind)=mean(edists);
    sedists(i1d,simind)=std(edists)./sqrt(length(edists)-1);
    dist2sim(i1d,simind)=Meandist(1,i1,1); %% just take undistorted distance to true
    set(gca,'Fontsize',18)
end; % id1/ invmethods

figure;hold on;
for i1=1:length(plotinv)
    h1=errorbar(i1,medists(i1),sedists(i1))
    set(h1,'Marker',mstr(i1));
    set(h1,'Color',colstr(i1));
    set(h1,'Linewidth',4);
    set(h1,'Markersize',20)
end; % for i1
axis([0 3 -Inf Inf])
set(gca,'XtickLabels',strvcat(' ',invmethods{plotinv(1)},invmethods{plotinv(2)},' '))
xlabel('Method')
ylabel('Distortion (mm)')
set(gca,'Fontsize',18)

figure;

se2dsim=std(Meandist(1,plotinv,:),[],3)./sqrt(Nseed-1);
hold on;

for f=1:length(plotinv)
    h1=plot(f,dist2sim(f),['*' colstr(f)]);
    set(h1,'Markersize',20)

end;
title([simpos ': Classic distance to true sim'])
ylabel('mm')

set(gca,'Xtick',1:length(plotinv))
set(gca,'Xticklabel',invmethods(plotinv))
set(gca,'Fontsize',18);
axis([0 numel(plotinv)+1 0 max(dist2sim)+2])




