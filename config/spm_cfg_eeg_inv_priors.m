function priors = spm_cfg_eeg_inv_priors
% Configuration file to set up priors for M/EEG source reconstruction
%__________________________________________________________________________
% Copyright (C) 2010-2017 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_cfg_eeg_inv_priors.m 7110 2017-06-15 11:45:05Z guillaume $


D = cfg_files;
D.tag = 'D';
D.name = 'M/EEG datasets';
D.filter = 'mat';
D.num = [1 Inf];
D.help = {'Select the M/EEG mat files.'};

val = cfg_entry;
val.tag = 'val';
val.name = 'Inversion index';
val.strtype = 'n';
val.help = {'Index of the cell in D.inv where the forward model can be found and the results will be stored.'};
val.val = {1};

priorname = cfg_entry;
priorname.tag = 'priorname';
priorname.name = 'Prefix for prior';
priorname.strtype = 's';
priorname.num = [1 Inf];
priorname.val = {'priorset1'};
priorname.help = {'Prefix for prior directory'};

sensorlevel = cfg_entry;
sensorlevel.tag = 'sensorlevel';
sensorlevel.name = 'RMS Sensor level noise';
sensorlevel.strtype = 'r';
sensorlevel.help = {'Assuming white sensor noise of same magnitude on all channels. Units fT for MEG. eg 10fT/rt(Hz) white noise for 100Hz bandwidth would be 100fT rms'};
sensorlevel.val = {100};


FWHMmm = cfg_entry;
FWHMmm.tag = 'FWHMmm';
FWHMmm.name = 'Patch smoothness';
FWHMmm.strtype = 'r';
FWHMmm.num = [1 1];
FWHMmm.val = {[4]};
FWHMmm.help = {'FWHM on cortex in mm. NB. Ignored for EBB and IID. Square windowed for COH '};


popular=cfg_menu;
popular.tag='popular';
popular.name='prior types';
popular.help={' Select from popular priors: IID (min norm), COH (LORETA) or EBB (beamforming)'};
popular.labels={'IID','COH','EBB','sparseEBB'};
popular.values={'IID','COH','EBB','sparseEBB'};
popular.val={'IID'};


popularpars=cfg_branch;
popularpars.tag='popularpars';
popularpars.name='Popular priors';
popularpars.help={' Select from popular priors and kernel smoothness (no smoothing used for IID)'};
popularpars.val={popular,FWHMmm};


npatches = cfg_entry;
npatches.tag = 'npatches';
npatches.name = 'Number of randomly selected patches';
npatches.strtype = 'i';
npatches.num = [1 1];
npatches.val = {[512]};
npatches.help = {'Number of randomly centred patches (priors) on each iteration'};

niter = cfg_entry;
niter.tag = 'niter';
niter.name = 'Number of iterations';
niter.strtype = 'i';
niter.num = [1 1];
niter.val = {[1]};
niter.help = {'Number of iterations'};


probfile = cfg_entry;
probfile.tag = 'probfile';
probfile.name = 'Prior to bias distribution of patches X';
probfile.strtype = 's';
%probfile.num = [1 1];
probfile.val = {''};
probfile.help = {'Name of patch probability file'};



randpatch = cfg_branch;
randpatch.tag = 'randpatch';
randpatch.name = 'Random Patches';
randpatch.help = {'Define random patches'};
randpatch.val  = {npatches,niter,FWHMmm,probfile};

fixedfile = cfg_files;
fixedfile.tag = 'fixedfile';
fixedfile.name = 'Patch definition file ';
fixedfile.filter = 'mat';
fixedfile.num = [1 1];
fixedfile.help = {'Select patch definition file (mat file with variable Ip: rows are iterations columns are patch indices '};
fixedfile.val={{''}};

fixedrows = cfg_entry;
fixedrows.tag = 'fixedrows';
fixedrows.name = 'Rows from fixed patch file';
fixedrows.strtype = 'i';
fixedrows.num = [1 2];
fixedrows.val = {[1 Inf]};
fixedrows.help = {'Select first and last row of patch sets to use from file'};


fixedpatch = cfg_branch;
fixedpatch.tag = 'fixedpatch';
fixedpatch.name = 'Fixed Patches';
fixedpatch.help = {'Define fixed patches'};
fixedpatch.val  = {fixedfile,fixedrows,FWHMmm};


sparsify = cfg_entry;
sparsify.tag = 'sparsify';
sparsify.name = 'Sparsify to N largest local maxima';
sparsify.strtype = 'i';
sparsify.num = [1 2];
sparsify.val = {[1 Inf]};
sparsify.help = {'Select first and last row of patch sets to use from file'};



space = cfg_menu;
space.tag = 'space';
space.name = 'Prior image space';
space.help = {'Space of the mask image.'};
space.labels = {'MNI', 'Native'};
space.values = {1, 0};
space.val = {1};

clean = cfg_menu;
clean.tag = 'clean';
clean.name = 'Start fresh';
clean.help = {'Start from fresh or keep existing priors'};
clean.labels = {'Fresh', 'Keep'};
clean.values = {1, 0};
clean.val = {1};


priorsmask  = cfg_files;
priorsmask.tag = 'priorsmask';
priorsmask.name = 'Priors file';
priorsmask.filter = '(.*\.gii$)|(.*\.mat$)|(.*\.nii(,\d+)?$)|(.*\.img(,\d+)?$)';
priorsmask.num = [0 1];
priorsmask.help = {'Select a mask or a mat file with priors or posterior.'};
priorsmask.val = {{''}};

funcpriors = cfg_branch;
funcpriors.tag = 'funcpriors';
funcpriors.name = 'Functionally defined priors';
funcpriors.help = {'Use posterior from previous inversion or restrict solutions to pre-specified VOIs'};
funcpriors.val  = {priorsmask, space,sparsify};


priortype = cfg_repeat;
priortype.tag = 'priortype';
priortype.name = 'Prior';
priortype.help = {'Specify the prior type to add'};
priortype.num  = [1 Inf];
priortype.values  = {randpatch,fixedpatch,funcpriors,popularpars};
priortype.val  = {randpatch};


%
% npost = cfg_entry;
% npost.tag = 'npost';
% npost.name = 'N posteriors';
% npost.strtype = 'n';
% npost.help = {'Nunber of posterior files to generate'};
% npost.val = {1};
%
% noisepost = cfg_entry;
% noisepost.tag = 'noisepost';
% noisepost.name = 'Percentage max noise';
% noisepost.strtype = 'r';
% noisepost.help = {'Add Gaussian noise with of a certain percentage of maximum power to posterior'};
% noisepost.val = {0};
%
%
% smoothpost = cfg_entry;
% smoothpost.tag = 'smoothpost';
% smoothpost.name = 'Number of smoothing iterations';
% smoothpost.strtype = 'r';
% smoothpost.help = {'Number of smoothing iterations (0 for no smoothing)'};
% smoothpost.val = {0};
%
% noisepost = cfg_entry;
% noisepost.tag = 'noisepost';
% noisepost.name = 'Percentage max noise';
% noisepost.strtype = 'r';
% noisepost.help = {'Add Gaussian noise with of a certain percentage of maximum power to posterior'};
% noisepost.val = {0};
%
% threshpost = cfg_entry;
% threshpost.tag = 'threshpost';
% threshpost.name = 'Percentage of maximum power at at which to threshold';
% threshpost.strtype = 'r';
% threshpost.help = {'Save posterior with full (0) or thresholded power distribution'};
% threshpost.val = {0};
%
% nmaxpost = cfg_entry;
% nmaxpost.tag = 'nmaxpost';
% nmaxpost.name = 'Number of maxima to keep as patch centers';
% nmaxpost.strtype = 'r';
% nmaxpost.help = {'Generate another file containing patch centres in Ip structure (0 make no file)'};
% nmaxpost.val = {0};
%
%
% postname = cfg_entry;
% postname.tag = 'postname';
% postname.name = 'Prefix for posterior';
% postname.strtype = 's';
% postname.num = [1 Inf];
% postname.val = {'postset1'};
% postname.help = {'Prefix for posterior directory'};
%


priors = cfg_exbranch;
priors.tag = 'priors';
priors.name = 'Inversion priors';
priors.val = {D, val,sensorlevel,priortype,priorname,clean};
priors.help = {'Set priors for source reconstruction'};
priors.prog = @add_priors;
priors.vout = @vout_priors;

%==========================================================================
function  out = add_priors(job)

val=job.val;
D = spm_eeg_load(job.D{val});

if ~isfield(D.inv{job.val}.inverse,'L')
    error('Need to set up spatial modes first');
else
    UL=D.inv{val}.inverse.L; %% dimension reduced lead field matrix
end


[a1,b1,c1]=fileparts(D.fname);
priordir=[D.path filesep job.priorname '_' b1];

if exist(priordir,'dir') ~= 7
    fprintf('Making directory for priors: %s\n',priordir);
    mkdir(priordir);
else
    mesh=D.inv{job.val}.mesh.tess_mni;
    [priorfiles] = spm_select('FPListRec',priordir,'.*\.mat$');
    
    count=size(priorfiles,1);
    fprintf('Found %d existing priors in this directory\n',count);
    
    if job.clean
        fprintf('Deleting existing priors\n');
        for f=1:count
            delete(deblank(priorfiles(f,:)));
        end
    else
        fprintf('Keeping existing priors\n');
    end
    
end

mesh=D.inv{val}.forward.mesh; %% assume this is in metres


Ns=size(mesh.vert,1);


% Set up sensor level prior
Qe{1}=eye(size(UL,1)).*job.sensorlevel;
Qp={}; %% source level covariance or svd components


% add in functional (from other modalities or experiment) hypotheses



for i=1:numel(job.priortype) %% move through different prior specifications
    allIp=[];
    Qp=[];
    
    % RANDOM PATCH CENTRES
    if isfield(job.priortype{i},'randpatch')
        
        base=job.priortype{i}.randpatch;
        
        Np = base.npatches;
        Npatchiter = base.niter;
        probfile=base.probfile;
        
        if isempty(probfile)
            prob=ones(Ns,1)./Ns;
        else
            
            probmesh=gifti(probfile);
            
            data=double(probmesh.cdata);
            M0.vertices=probmesh.vertices;
            M0.faces=probmesh.faces;
            
            if size(probmesh.vertices,1)==size(mesh.vert,1)
            %    fprintf('NB smoothing prior with 20 iterations\n');
                %prob=spm_mesh_smooth(M0,data',20);
                prob=data;
            else
                error('Gifti mesh size does not match');
            end
        end
        if sum(prob)>1.01
            warning('Rescaling bias file to unit sum');
        end
        
        prob=prob./sum(prob);
        Ip=[];
        
        
        fprintf('Using %d iterations of %d random patches\n',Npatchiter,Np);
        allIp=sparse(Npatchiter,Np);
        
        
        %rng('shuffle');
        
        dum=sparse(1,Ns);
        figure;
        hold on;
        
        for j=1:Npatchiter
            
            for k=1:Np
                pval=max(prob)*(k-1)./Np;
                sind=find(prob>=pval);
                randind=randperm(length(sind));
                allIp(j,k)=sind(randind(1));
                dum(allIp(j,k))=dum(allIp(j,k))+1;
            end
            %end; %% IF
            plot(1:length(prob),dum,'.');  drawnow;
        end % for j
        fprintf('Sampled %3.2f percent of cortex (%d vertices) \n',100*length(find(dum>0))./length(dum),length(find(dum>0)));
        
    end % random patch
    
    % FIXED PATCH CENTRES FROM FILE
    if isfield(job.priortype{i},'fixedpatch')
        
        base=job.priortype{i}.fixedpatch;
        try
            dum=load(char(base.fixedfile));
        catch
            error('failed to load fixedpatch file');
        end
        if ~isfield(dum,'Ip')
            error('Need to have patch indices in structure Ip');
        end
        allIp=dum.Ip;
        rows=base.fixedrows;
        lastrow=min(rows(2),size(allIp,1));
        rowind=[rows(1):lastrow];
        allIp=allIp(rowind,:);
    end % fixed patch
    
    %%%%%%%%%%%%%%% FIXED OR RANDOM PATCHES
    if isfield(job.priortype{i},'fixedpatch') || isfield(job.priortype{i},'randpatch')
        
        % Set up patches on cortical surface determined by indices of allIp
        % one prior files per row of Ip
        
        [Qp,Qe,allpriornames]=spm_eeg_invert_setuppatches(allIp,mesh,base,priordir,Qe,UL);
        
        
        
    end % if fixedpatch or randpatch
    
    
    % figure;
    %
    % if job.nmaxpost>0,
    %     postIpname=[postdir filesep sprintf('post_Ip.mat')];
    %     Ip=zeros(job.npost,job.nmaxpost);
    % end;
    %
    % M0=[];M0.vertices=mesh.vert; M0.faces=mesh.face;
    % for f=1:job.npost,
    %
    %     subplot(job.npost,1,f);
    %     plot(sourcevar,'k');hold on;
    %     postname=[postdir filesep sprintf('post%d.mat',f)];
    %     fprintf('adding %f percent noise to posterior\n',job.noisepost);
    %     nsourcevar=sourcevar+max(sourcevar)*job.noisepost*randn(size(sourcevar))./100;
    %     fprintf('Smoothing by %d iterations\n',job.smoothpost);
    %     nsourcevar = spm_mesh_smooth(M0, nsourcevar',job.smoothpost);
    %
    %     fprintf('Thresholding at %f percent of max power\n',job.threshpost);
    %     thresh=max(nsourcevar)*job.threshpost./100;
    %     plot([0,length(sourcevar)],[thresh thresh],':g');
    %     nsourcevar(nsourcevar<thresh)=0;
    %     plotind=find(nsourcevar);
    %     plot(plotind,nsourcevar(plotind),'ro');
    %     legend('orig','thresh','thresh noisy');
    %     nsourcevar=sparse(nsourcevar);
    %     Qp{1}=diag(nsourcevar);
    %     fprintf('SAving %s\n',postname);
    %     F=-Inf;
    %     save(postname,'Qp','Qe','UL','F','-v7.3');
    %     %%
    %     fprintf('\n Getting %d local maxima ',job.nmaxpost);
    %     Ip1 = spm_mesh_get_lm(M0,nsourcevar); %% get local maxima
    %
    %     [vals,ind]=sort(nsourcevar(Ip1),'descend');
    %     Ip(f,:)=Ip1(ind(1:job.nmaxpost));
    %
    %
    %
    % end; % for f
    % if job.nmaxpost>0,
    %     fprintf('Saving %s\n',postIpname);
    %     save(postIpname,'Ip'); % =[postdir filesep sprintf('post_Ip.mat')];
    %     %Ip=zeros(job.npost,job.nmaxpost);
    % end;
    
    
    
    
    %%%%%%%%%%% FUNCTIONALLY DEFINED PRIORS
    
    if isfield(job.priortype{i},'funcpriors')
        base=job.priortype{i}.funcpriors;
        P = char(base.priorsmask);
        priors{i}.priortype='funcpriors';
        priors{i}.priorfname=P;
        
        fprintf('Using priors from file %s\n',P);
        
        P=priors{i}.priorfname;
        [p,f,e] = fileparts(P);
        switch lower(e)
            case '.mat'
                a=load(P);
                if isfield(a,'Qp')
                    
                    disp('Replacing sensor and source posterior distribution');
                    for f=1:numel(a.Qp)
                        Qp{end+1}=a.Qp{f};
                    end
                    
                    % now need to compress the original Qe (say 275 channels) into reduced channel form
                    U=D.inv{val}.inverse.U{1}; %% get spatial modes definition
                    disp('NOT SURE IF WE NEED TO RESCALE FOR CHAN REDUCTION');
                    scaleU=size(a.Qe,1)./size(U,1); % white noise power pwer channel reduce with number of channels
                    
                    Qe{1}=U*a.Qe*U'./scaleU;
                    fprintf('\n Reducing white noise variance from %3.2f fT^2 on each of %d channels to %3.2f fT^2 on each of %d channels\n',...
                        mean(diag(a.Qe)),size(a.Qe,1),mean(diag(Qe{1})),size(U,1));
                    
                    
                else
                    disp('loading just source level posterior');
                    Qp{end+1}=a.pQ;
                end
                
                
                
            case '.gii'
                g = gifti(P);
                
                for i=1:size(g.cdata,2)
                    Qp{end+1} = double(g.cdata(:,i));
                end
                
                
            case {'.img', '.nii'}
                S.D = D;
                S.fmri = P;
                S.space = job.isstandard.custom.priors.space;
                D = spm_eeg_inv_fmripriors(S);
                inverse.fmri = D.inv{D.val}.inverse.fmri;
                a=load(inverse.fmri.priors);
                
                Qp{end+1} = a.pQ;
            otherwise
                error('Unknown file type.');
        end % switch
        
        idnum=randi(1e6);
        
        priorfname=[priordir filesep sprintf('prior%d.mat',idnum)];
        fprintf('Saving %s\n',priorfname);
        F=[]; % no associated free energy value
        save(priorfname,'Qp','Qe','UL','F');
    end % func priors
    
    % POPULAR PRIORS
    
    if isfield(job.priortype{i},'popularpars')
        base=job.priortype{i}.popularpars;
        
        fprintf('Adding prior covariance matrix for %s algorithm\n',base.popular);
        
        if strcmp(base.popular,'sparseEBB')
            
            % calculate power on cortical surface
            % using beamformer assumptions
            AYYA=D.inv{val}.inverse.AYYA;
            
            InvCov = spm_inv(AYYA);
            
            allsource = sparse(Ns,1);
            Sourcepower = sparse(Ns,1);
            
            for bk = 1:Ns
                normpower = 1/(UL(:,bk)'*UL(:,bk));
                Sourcepower(bk) = 1/(UL(:,bk)'*InvCov*UL(:,bk));
                % normalise power by unit noise through lead fields
                allsource(bk) = Sourcepower(bk)./normpower;
            end
            
            % now get local maxima on mesh
            M.vertices=mesh.vert;
            M.faces=mesh.face;
            
            Ip = spm_mesh_get_lm(M,allsource); %% get local maxima
            figure;
            plot(allsource);
            legend('sparse EBB');
            hold on;
            maxBFpatch=100;
            fprintf('Limiting to a max of %d peaks\n',maxBFpatch);
            
            [vals,ind]=sort(allsource(Ip),'descend');
            Ip=Ip(ind(1:maxBFpatch));
            plot(Ip,allsource(Ip),'ro');
            
            [Qp,Qe,allpriornames]=spm_eeg_invert_setuppatches(Ip,mesh,base,priordir,Qe{1},UL);
            
            
        end % sparse EBB
        
        if strcmp(base.popular,'EBB')
            
            % calculate power on cortical surface
            % using beamformer assumptions
            AYYA=D.inv{val}.inverse.AYYA;
            
            InvCov = spm_inv(AYYA);
            
            allsource = sparse(Ns,1);
            Sourcepower = sparse(Ns,1);
            
            for bk = 1:Ns
                normpower = 1/(UL(:,bk)'*UL(:,bk));
                Sourcepower(bk) = 1/(UL(:,bk)'*InvCov*UL(:,bk));
                % normalise power by unit noise through lead fields
                allsource(bk) = Sourcepower(bk)./normpower;
            end
            
            figure;
            plot(allsource);
            legend('EBB');
            hold on;
            % now get local maxima on mesh
            M.vertices=mesh.vert;
            M.faces=mesh.face;
            warning('Smoothing not used');
            Qp{1}=diag(allsource);
            
            
            
            
            idnum=randi(1e6);
            priorfname=[priordir filesep sprintf('prior%d.mat',idnum)];
            F=[];
            save(priorfname,'Qp','Qe','UL','F');
            
        end % EBB
        
        if strcmp(base.popular,'COH')
            
            % create minimum norm prior
            %------------------------------------------------------------------
            Qp{1}  = speye(Ns,Ns);
            warning('SQUARE SMOOTHING KERNEL OF FWHM MM');
            
            M=[];M.faces=mesh.face;M.vertices=mesh.vert;
            dist=spm_mesh_distmtx(M,1);
            ind2=find(dist>0);
            ind=find(dist<base.FWHMmm./1000);
            ind=intersect(ind,ind2);
            Qp{1}(ind)=1;
            
            disp('Set up COH prior');
            
            idnum=randi(1e6);
            priorfname=[priordir filesep sprintf('prior%d.mat',idnum)];
            fprintf('Saving %s\n',priorfname);
            F=[]; % no associated free energy value
            
            save(priorfname,'Qp','Qe','UL','F');
        end % COH
        
        if strcmp(base.popular,'IID')
            
            % create minimum norm prior
            %------------------------------------------------------------------
            Qp{1}  = speye(Ns,Ns);
            
            disp('Set up IID prior');
            
            idnum=randi(1e6);
            
            priorfname=[priordir filesep sprintf('prior%d.mat',idnum)];
            fprintf('Saving %s\n',priorfname);
            F=[]; % no associated free energy value
            save(priorfname,'Qp','Qe','UL','F');
        end % IID
        
        
    end % popular priors
    
end % for i (looping over priors)



% NOW JUST WANT TO PLOT OUT THE CURRENT DENSITY PRIORS FOR A RANDOM SET OF PATCHES
% THIS SHOULD ALL GO IN A PREVIEW FUNCTION

[LCpL,Q,sumLCpL,QE,Cy,M,Cp,Cq,Lq]=spm_eeg_assemble_priors(UL,Qp,Qe);


figure;
mesh=D.inv{job.val}.mesh.tess_mni;
M0.vertices=mesh.vert;
M0.faces=mesh.face;
sourcepwr=diag(Cp);

spm_mip(sourcepwr,mesh.vert',6);
title('Prior');
colorbar;


CHECKSMOOTH=0;
if CHECKSMOOTH
    figure;
    r=randperm(length(Qp));
    
    sigma2=[];
    for f=1:min(length(r),4)
        if isfield(Qp{r(f)},'q')
            q=Qp{r(f)}.q;
            [val,ind]=max(q);
            
            order=1; %1 distance in mm, 0 distance in edges
            dist = spm_mesh_geodesic(M,ind-1,order);
            distind=intersect(find(dist<0.020),find(q>0)); % less than 2cm away
            
            y=spm_vec(log(q(distind)'));
            x=spm_vec(-dist(distind).^2);
            B = regress(y,[x ones(size(x,1),1)] );
            sigma2(f)=1./(2*B(1)); %% sd squared of gaussian
            plot(dist(distind).*1000,q(distind),'r.',dist(distind).*1000,max(q).*exp(-(dist(distind).^2)/(2*sigma2(f))),'go'); %% plot in mm
        end
        
    end
    
    if ~isempty(sigma2)
        title(sprintf('Prior Profile of prior amplitude across cortex'));
        meanfwhm=mean(sqrt(sigma2))*2.355*1000; %% in mm
        stdfwhm=std(sqrt(sigma2))*2.355*1000;
        legend(sprintf('Mean FWHM=%3.2f mm,std FWHM= %3.2fmm',meanfwhm,stdfwhm));
        drawnow;
    end
    
end %% if CHECKSMOOTH

out.D = job.D;

%==========================================================================
function dep = vout_priors(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) after imaging source reconstruction';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});


% function [q,dist,useind]= gauss_patch(M,i,FWHM,q)
%
% order=1;
%
%
% sigma=FWHM./2.355;
% sigma2=sigma^2;
%
% d=spm_mesh_geodesic(M,i-1,order);
% useind=find(d<FWHM*2);
% dist=d(useind);
% q(useind)=exp(-(dist.^2)/(2*sigma2));
