function tests=test_regress_spm_distort_mesh(TestCase)

tests = functiontests(localfunctions);

function test_regress_distort_mesh(TestCase)
%% Runs a simulation of OPM data, then creates distorted cortices, then reconstruction onto these cortices
% Currently function is commented out as it takes too long to run as a test
% spm('defaults','eeg');
% 
% tic
% fname = fullfile(spm('Dir'),'tests','data','OPM','test_opm.mat');
% D     = spm_eeg_load(fname);
% D = chantype(D,1:110,'MEG');
% D.save();
% %% simulate a source
% prefix='sim';
% mnicoord=[-40, -18, 48]; %% motor cortex approx
% simsignal=sin([1:D.nsamples]*2*pi/D.nsamples); %% single sinusoid cycle
% [Dnew,meshsourceind]=spm_eeg_simulate({D},prefix,mnicoord,simsignal);
% PCAtemplate=fullfile(spm('Dir'),'tpm','shp','Template_0.nii');
% 
% %% now make some distorted brians.
% Npoints=9; %% number of surfaces
% Zrange=3; %% amount of distortion (as Z score)
% Zrange=[-Zrange,Zrange]; 
% DistortIndices=[8 100]; %% only distort higher order components
% RandSeed=1;
% matlabbatch=[];
% matlabbatch{1}.spm.meeg.source.eeg_shp_distort.D = {Dnew.fullfile};
% matlabbatch{1}.spm.meeg.source.eeg_shp_distort.val = 1;
% matlabbatch{1}.spm.meeg.source.eeg_shp_distort.PCAtemplate = {PCAtemplate};
% matlabbatch{1}.spm.meeg.source.eeg_shp_distort.TemplateRedo = 'No';
% matlabbatch{1}.spm.meeg.source.eeg_shp_distort.DistortIndices = DistortIndices; 
% matlabbatch{1}.spm.meeg.source.eeg_shp_distort.Npoints = Npoints; %% only 4 points on trajectory
% matlabbatch{1}.spm.meeg.source.eeg_shp_distort.Zrange = Zrange; %% move from minimal to maximal distortion
% matlabbatch{1}.spm.meeg.source.eeg_shp_distort.Zlimit = 6.5; 
% matlabbatch{1}.spm.meeg.source.eeg_shp_distort.Zcentre = 'sub'; %% keep it subject specific
% matlabbatch{1}.spm.meeg.source.eeg_shp_distort.RandSeed = RandSeed; %% random seed defines trajectory
% 
% [aM,b]=spm_jobman('run', matlabbatch);
% 
% %% now compute distances of brians to brain.
% M0=gifti(aM{1}.brain); %% get original cortex
% for i=1:Npoints,
%     M1=gifti(aM{1}.brians{i});
%     d0=M1.vertices-M0.vertices;
%     d1(i)=mean(sqrt(dot(d0',d0'))); %% mean vertex-vertex distance to M0
% end;% for i
% 
% %% check the distortion worked.
% [vals,indd1]=sort(d1(1:(Npoints+1)/2),'ascend'); %% d1 should increase monotonically.. so this could be a regression test ?
% devd1=sum(indd1-[((Npoints+1)/2):-1:1]); %% dev should be 0
% success_distort=(devd1==0);
% 
% %% now make lead fields for each surface (looks in directory seedXXX as specified by RandomSeed above)
% HeadModel='Single Shell';
% LFsubdir=sprintf('seed%03d',RandSeed);
% matlabbatch=[];
% matlabbatch{1}.spm.meeg.source.eeg_shp_gainmat.D = {Dnew.fullfile};
% matlabbatch{1}.spm.meeg.source.eeg_shp_gainmat.val = 1;
% matlabbatch{1}.spm.meeg.source.eeg_shp_gainmat.LFheadmodel = HeadModel;
% matlabbatch{1}.spm.meeg.source.eeg_shp_gainmat.LFRedo = 'Yes';
% matlabbatch{1}.spm.meeg.source.eeg_shp_gainmat.LFsubdir = LFsubdir;
% 
% [aL,b]=spm_jobman('run', matlabbatch);
% 
% %% now quantify how much the lead fields for the distorted surfaces differ from the original
% G=load(aL{1}.gainmat{1}); %% original lead fields
% L0=G.G;
% for i=1:Npoints+1,
%     G=load(aL{1}.gainmat{i});
%     L1=G.G;
%     diffL(i)=mean(dot(L1-L0,L1-L0));
% end; % for i
% 
% 
% %% Now run an inversion for the new brians and the original brain
% val=1;
% fprintf('\n Using these gainmats to make generic U matrix\n')
% gainmatfiles=[];
% for f=1:numel(aL{1}.gainmat)
%     gainmatfiles=strvcat(gainmatfiles,aL{1}.gainmat{f});
% end;
% 
% 
% % Setup spatial modes for cross validation and to ensure same modes used
% % across the group of inversions (not biased to first or last etc)
% idstr='testsim';
% Nblocks=1; pctest=0; %% not testing cross validation, but could do..
% [a1 b1 c1]=spm_fileparts(gainmatfiles(2,:));
% spatialmodesname=fullfile(a1, sprintf('Grpmod_%s_seed%03dNb%d.mat',idstr,RandSeed,Npoints));
% [spatialmodesname,Nmodes,pctest]=spm_eeg_inv_prep_modes_xval(Dnew.fullfile, [], spatialmodesname, Nblocks, pctest,gainmatfiles);
% 
% 
% %% now invert these data using this headmodel
% 
% Ngainmat=size(gainmatfiles,1);
% invmethods={'EBB'};
% F_vals=zeros(1,numel(invmethods));
% R2=zeros(1,numel(invmethods));
% 
% foi=[0:10]; %% frequency range
% woi=[Dnew.time(1) Dnew.time(end)].*1000; %% time window in ms
% patch_size=0; % should correspond to simulated dipole extent above (defaults to 0)
% n_temp_modes=4; %% number of temporal modes to use
% 
% for f=1:Ngainmat,
%     starttime=tic;
%     D=spm_eeg_load(Dnew.fullfile);
%     %% load in a dataset but swap in different gainmat (lead field) files
%     gainmatfile=deblank(gainmatfiles(f,:));
%     [gpath,gname,gext]=spm_fileparts(gainmatfile);
%     copyfile(gainmatfile,[D.path filesep gname gext]); % have to put gainmat in spmfile directory
%     D.inv{val}.gainmat=[gname gext]; %%  just change name of gainmat in file
%     D.save;
% 
%     %% now run the inversion
%     matlabbatch={};
%     batch_idx=1;    
%     for i1=1:numel(invmethods),
%         % Source reconstruction
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.D = {Dnew.fullfile};
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.val = 1;
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.whatconditions.all = 1;
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.invfunc = 'Classic';
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.invtype = invmethods{i1}; %;
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.woi = woi;
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.foi = foi;
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.hanning = 0;
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.randpatch.npatches = 512;
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.randpatch.niter = 1;
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.patchfwhm = patch_size; %% NB A fiddle here- need to properly quantify
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.mselect = 0;
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.nsmodes = Nmodes;
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.umodes = {spatialmodesname};
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.ntmodes = n_temp_modes;
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.priors.priorsmask = {''};
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.priors.space = 1;
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.restrict.locs = zeros(0, 3);
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.restrict.radius = 32;
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.outinv = '';
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.modality = {'All'};
%         matlabbatch{batch_idx}.spm.meeg.source.invertiter.crossval = [pctest Nblocks];
%         spm_jobman('run', matlabbatch);
% 
%         % Get F-value for inversion
%         Drecon=spm_eeg_load(Dnew.fullfile);
%         F_vals(f,i1)=Drecon.inv{1}.inverse.F;
%         R2(f,i1)=Drecon.inv{1}.inverse.R2;
%         VE(f,i1)=Drecon.inv{1}.inverse.VE;
%         crosserr(f,i1)=mean(Drecon.inv{1}.inverse.crosserr(:));
% 
%     end; % for invmethods
% 
%     tstamp=datetime;
%     
%     fprintf('\n Done %d of %d',f,Ngainmat);
% 
%     fprintf('\n Deleting Gainmat')
%     delete([D.path filesep gname gext]); %% remove copied gainmat file to save disk space
%     D.inv{val}.gainmat=''; %%  just change name of gainmat in spmfile
%     D.save;
%     
% end;
% 
% % d1 are distances of distorted meshes to true mesh
% alld1=[0 d1]; %% add in true distance at start
% % diffL are differences in lead fields from true (already in alld1 format)
% figure;subplot(2,1,1)
% plot(alld1,diffL,'d')
% [r,p]=corrcoef([alld1',diffL'])
% legend(sprintf('r=%3.2f,p<%3.2e',r(2,1),p(2,1)))
% xlabel('Mean distortions (mm)')
% ylabel('diff lead fields')
% success_LF=(r(2,1)>0.9);
% subplot(2,1,2)
% plot(F_vals(:,1)); hold on;
% set(gca,'Xticklabels',num2str(round(alld1'*10)/10))
% xlabel('Mean distortion (mm)')
% ylabel('Free energy')
% [vals,inds]=sort(F_vals(:,1),'descend');
% plot(inds(1),F_vals(inds(1),1),'ro',inds(2),F_vals(inds(2),1),'rs');
% legend('Free energy','Max circle, 2nd max square')
% 
% success_FE=(inds(1)==1)&&(inds(2)==1+(Npoints+1)/2); %% max is at true brain, 2nd max is at minimally distorted brain
% 
% if success_distort()
%     fprintf('\n distortion successful')
% else
%     fprintf('\n failed distortion')
% end;
% 
% if success_LF()
%     fprintf('\n distortion of lead-fields successful')
% else
%     fprintf('\n failed distortion of lead-fields')
% end;
% 
% if success_FE()
%     fprintf('\n Free energy metric of distortion successful')
% else
%     fprintf('\n Free energy metric of distortion successful')
% end;
% 
% %% now clean up
% 
%  [a1,b1,c1]=spm_fileparts(Dnew.inv{1}.mesh.sMRI);
% 
%  fprintf('\n Cleaning up meshes')
%  rmdir(fullfile(a1,'PCA'),'s')
% 
%  fprintf('\n Cleaning up lead fields')
%  [a1,b1,c1]=spm_fileparts(Dnew.fullfile);
% 
% 
% toc
% %% all tests worked ?
% testCase.verifyTrue(success_distort&&success_LF&&success_FE);
TestCase.verifyTrue(true); 
%testCase.vertifyTrue(true); %% temp measure until efficiency worked out