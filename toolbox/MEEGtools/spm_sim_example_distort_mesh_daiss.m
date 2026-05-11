clear all;
close all;

spm('defaults','eeg');




%% after that =0 just loads in results.

fname = fullfile(spm('Dir'),'tests','data','meg_face','cdbespm8_SPM_CTF_MEG_example_faces1_3D.mat');


if ~exist(fname),
    error('Need to download SPM test data')
end;

[datadir,b1,c1]=fileparts(fname);

D     = spm_eeg_load(fname);



%% simulate a source
prefix='sim';



HeadModel='Single Shell'; % forward model to use


%% now make some distorted brians.
Npoints=9; % 17; %% number of surfaces
Zrange=3; %% amount of distortion (as Z score)
Nseed=8; %% number of random seeds to use
Zrange=[-Zrange,Zrange];
DistortIndices=[8 100]; %% only distort higher order components




simstr={'1motor'}; %% simulation scenario


Jpk=[];
Fieldmap=[];
Lf=[];




tic

%% simulate a source
prefix='sim';

RUNSIM=1; %% run all the simulations (1) or just load in data (0)
tic



simpos='1motor'
switch (simpos),
    case '1motor',
        mnicoord=[-40, -18, 48]; %% motor cortex approx
        freq=30; %% frequency of sinusoid
        scale=[1]; % how much to scale simulated source by

        %% Can add other scenarios here
    otherwise
        error('case not found')

end; % case simpos

foi=[0 freq*2];
invmethods={'eLOR','sLOR','LCMV'};


if RUNSIM,
    mult_simsignal=sin(freq*[1:D.nsamples]*2*pi/D.nsamples).*scale; %% single sinusoid cycle
    [Dnew,meshsourceind]=spm_eeg_simulate({D},prefix,mnicoord,mult_simsignal);
    PCAtemplate=fullfile(spm('Dir'),'tpm','shp','Template_0.nii');

    %% now make some distorted brians.
    
    
    
    allBFnames=[];
    for RandSeed=1:Nseed,
        matlabbatch=[];
        matlabbatch{1}.spm.meeg.source.eeg_shp_distort.D = {Dnew.fullfile};
        matlabbatch{1}.spm.meeg.source.eeg_shp_distort.val = 1;
        matlabbatch{1}.spm.meeg.source.eeg_shp_distort.PCAtemplate = {PCAtemplate};
        matlabbatch{1}.spm.meeg.source.eeg_shp_distort.TemplateRedo = 'No';
        matlabbatch{1}.spm.meeg.source.eeg_shp_distort.DistortIndices = DistortIndices;
        matlabbatch{1}.spm.meeg.source.eeg_shp_distort.Npoints = Npoints; %% N points on trajectory
        matlabbatch{1}.spm.meeg.source.eeg_shp_distort.Zrange = Zrange; %% move from minimal to maximal distortion
        matlabbatch{1}.spm.meeg.source.eeg_shp_distort.Zlimit = 6.5;
        matlabbatch{1}.spm.meeg.source.eeg_shp_distort.Zcentre = 'sub'; %% keep it subject specific
        matlabbatch{1}.spm.meeg.source.eeg_shp_distort.RandSeed = RandSeed; %% random seed defines trajectory

        [aM,b]=spm_jobman('run', matlabbatch);

        %% now compute distances of brians to brain.
        M0=gifti(aM{1}.brain); %% get original cortex
        for i=1:Npoints,
            M1=gifti(aM{1}.brians{i});
            d0=M1.vertices-M0.vertices;
            d1(i)=mean(sqrt(dot(d0',d0'))); %% mean vertex-vertex distance to M0
        end;% for i


        %% now run DAISS code for each surface.


        nsHeadModel=replace(HeadModel,' ','_')
        LFsubdir=sprintf('seed%03d',RandSeed);
        val=1;





        for i=1: Npoints+1,
            if i==1,
                cmesh=aM{1}.brain; %% true anatomy
            else
                cmesh=aM{1}.brians{i-1}; % surrogate anatomy
            end;
            M=Dnew.inv{1}.datareg.toMNI;
            Dnew.inv{1}.forward.toMNI=M;

            Mc=gifti(cmesh);
            %%
            warning('No deformation used')
            mnivert= spm_eeg_inv_transform_points(inv(M), Mc.vertices);
            %% need to distort in mni coords for daiss
            Dnew.inv{val}.mesh.tess_mni.vert=mnivert.*1000;
            Dnew.save;

            for invind=1:numel(invmethods),
                invmethod=invmethods{invind};
                for oind=1:2, %% orientation
                    if oind==1,
                        orientstr='original';
                    else
                        orientstr='unoriented';
                    end;

                    matlabbatch=jimmydaiss_job_batch_bf(Dnew,orientstr,foi,invmethod,datadir)
                    
                    [a,b]=spm_jobman('run', matlabbatch);
                    BFname=a{1}.BF;
                    BF=load(BFname{1});
                    [a1,b1,c1]=fileparts(BFname{1});
                    saveBFname=fullfile(a1,sprintf('Oct_BFpt%d_%s_orient%d_seed%03d_%s_%d.mat',i,invmethods{invind},oind,RandSeed,nsHeadModel,round(scale*1000)));
                    save(saveBFname,'BF');
                    allBFnames=strvcat(allBFnames,saveBFname);

                    close all;
                end; % oind
            end; % invmethod



        end;
    end; % for randseed
    elapsed=toc;
    save([datadir filesep 'simwkspce.mat']);
else % if RUNSIM
    load([datadir filesep 'simwkspce.mat']);
end; % if RUNSIM

[outdir,b1,c1]=fileparts(allBFnames(1,:));




close all;
for invind=1:numel(invmethods),
    invmethod=invmethods{invind};
    for oind=1:2, %% orientation
        if oind==1,
            orientstr='original';
        else
            orientstr='unoriented';
        end;
        maxpwrvals=zeros(Nseed,Npoints+1);
        for RandSeed=1:Nseed,
            for i=1:Npoints+1,
                loadBFname=fullfile(a1,sprintf('Oct_BFpt%d_%s_orient%d_seed%03d_%s_%d.mat',i,invmethods{invind},oind,RandSeed,nsHeadModel,round(scale*1000)));
                %loadBFname=fullfile(outdir,sprintf('Oct_BFpt%d_%s_orient%d_seed%03d_%s.mat',i,invmethods{invind},oind,RandSeed,nsHeadModel));
                load(loadBFname,'BF');
                maxpwrvals(RandSeed,i)=max(BF.output.image.val);
            end; % for i
            RandSeed
        end; % RandSeed
        [shortPwrvals,shortdist]=simplify_dist_metric(maxpwrvals,[0 d1]);
        allpwrvals(oind,invind,:,:)=maxpwrvals;
        [dum,ind]=max(shortPwrvals');
        meandist(oind,invind)=mean(shortdist(ind));
        sddist(oind,invind)=std(shortdist(ind))./sqrt(Nseed-1);
        figure;
        plot(shortdist,shortPwrvals);
        title([orientstr ' ' invmethods{invind}]);
        hold on;
        plot(shortdist(ind),dum,'x');
    end;% orientation
end; % invmethod
figure; hold on;
mstr='sd'
legstr='';
for oind=1:2,
    xvals=[1:numel(invmethods)+(oind-1)/2];
    h=errorbar(xvals,meandist(oind,:),sddist(oind,:),mstr(oind))
    if oind==1,
        orientstr='original';
    else
        orientstr='unoriented';
    end;
    legstr=strvcat(legstr,orientstr)
    set(h,'LineWidth',4);
    set(h,'MarkerSize',10)
end;
legend(legstr)
set(gca,'Xtick',1:numel(invmethods));
set(gca,'XtickLabel',invmethods)
ylabel('Distortion (mm)')
xlabel('Algorithm')
title(['Based on peak projected power'])
set(gca,'Fontsize',18)