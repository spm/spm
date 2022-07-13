function [allsmoothnames] = spm_eeg_smoothmesh_mm(meshname,allfwhm,redo)
%%function [allsmoothnames] = spm_eeg_smoothmesh_mm(meshname,allfwhm,redo)
%% produces a smoothing kernel for mesh in mm - each colum QG(:,j)
%% if redo==1 and file already exists then redo

% $Id: spm_eeg_invert_setuppatches.m 6498 2015-07-15 19:13:31Z gareth $



[a1,b1,c1]=fileparts(meshname);
M=gifti(meshname);

Ns=size(M.vertices,1);
if nargin<3,
    redo=[];
end;

if isempty(redo)
    redo=0;
end;

A=spm_mesh_area(M,1);

allsmoothnames=[];
for k=1:length(allfwhm),
    fwhm=allfwhm(k);
    smoothmeshname=[a1 filesep sprintf('FWHM%3.2f_',fwhm) b1 '.mat']
    allsmoothnames=strvcat(allsmoothnames,smoothmeshname);
    if exist(smoothmeshname) && ~redo,
        fprintf('\nFound smoothfile %s, not recomputing\n',smoothmeshname);
        continue;
    end;
    vertspace=mean(sqrt(A));
    fprintf('\n FWHM of %3.2f is approx %3.2f times vertex spacing\n',fwhm, fwhm/vertspace);
    if (fwhm>vertspace*100) || (fwhm<vertspace/100),
        error('Mismatch between FWHM and mesh units');
    end;
    
    
    QG=zeros(Ns,Ns);
    sigma=fwhm./2.355;
    sigma2=sigma^2;
    
    order=1;
    
    fprintf('\n Smoothing with fwhm %3.2f (%d of %d) this make take some time (~20mins per FWHM with parallel toolbox)\n',fwhm,k,length(allfwhm))
    %parfor j = 1:Ns,
    for j = 1:Ns,
        % Patch locations determined by Ip
        %--------------------------------------------------------------
        q=zeros(1,Ns);
        d = spm_mesh_geodesic(M,j-1,order);
        
        useind=find(d<=fwhm); %% maybe only extend just over 2 sigma in any one direction (i.e. cover 95% interval)
        q(useind)=exp(-(d(useind).^2)/(2*sigma2));
        q=q.*(q>exp(-8));
        q=q./sum(q);
        QG(:,j)=q;
        
    end; % for j
    
    
    
    fprintf('Saving %s\n',smoothmeshname);
    
    
    QG=sparse(QG);
    save(smoothmeshname,'QG','M','-v7.3');
end; % for k





