function [Qp,Qe,allpriornames] = spm_eeg_invert_setuppatches(allIp,mesh,base,priordir,Qe,UL)
% Set prior files for source inversion
% FORMAT [Qp,Qe,allpriornames] = spm_eeg_invert_setuppatches(allIp,mesh,base,priordir,Qe,UL)
% Each file contains  number of smooth patches on cortical surface a
% allIp    - each row denotes a different prior file
%            each column denotes the index of an impulse on the cortical surface
% mesh     - cortical surface mesh (in metres)
% base.nAm (optional)    - magnitude of the impulse.
%                          There should be one value per column of Ip
% base.smooth (optional) - FWHM smoothness of the impulse on cortical surface (in mm)
% priordir - Directory in which the new priorfiles will be saved
% Qe       - sensor level covariance
% UL       - reduced lead field (only used to make a complete prior file)
%
% Qp  - prior source covariances from prior created in last row of allIp
% Qe  - prior sensor covariances
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_eeg_invert_setuppatches.m 7118 2017-06-20 10:33:27Z guillaume $


Npatchiter=size(allIp,1);
Np=size(allIp,2);%% number of patches per iteration
smoothm=base.FWHMmm./1000;
disp(sprintf('Using %d iterations of %d fixed patches',Npatchiter,Np));


M.vertices=mesh.vert;
M.faces=mesh.face;

Ns=size(M.vertices,1);

if isfield(base,'nAm')
    nAm=base.nAm;
else
    nAm=ones(Np,1).*10;
    fprintf('\nNo magnitudes defined, setting to %3.2fnAm\n',nAm(1));
end;

if isfield(base,'FWHMmm')
    smoothm=base.FWHMmm./1000;
else
    smoothm=5;
    fprintf('\nSmoothness of all patches is set to %3.2fmm\n',smoothm(1)*1000);
end;
if length(smoothm)==1,
    smoothm=ones(Np,1).*smoothm; %% turn it into a vector
end;


[priorfiles] = spm_select('FPListRec',priordir,'.*\.mat$');
priorcount=size(priorfiles,1);
allpriornames=[];

for k=1:Npatchiter,
    Ip=allIp(k,:);
    Qp={};
    
    for j = 1:Np
        % Patch locations determined by Ip
        %--------------------------------------------------------------
        
        q0=sparse(1,Ns).*0;
        
        [q,dist,useind,areapervertex]               = gauss_patch(M,Ip(j),smoothm(j),q0); %% priors.smooth should be in metres like the mesh
        
        q=q./sum(q); %patch adds up to unity
        q=q.*nAm(j);
        
        Qp{j}.q   = q';
        if length(dist)>1, %% more than 1 vertex in fwhm
            areamm2=pi*(max(dist*1000)/2).^2;
        else
            areamm2=areapervertex*1e6; % take approx area per vertex in this region
        end;
        
        mompervertex=mean(q(useind));
        peakmom=max(q(useind));
        mompermm2(j)=full(mompervertex*length(dist)/(areamm2)); %% nAm/mm2
        peakmompermm2(j)=full(peakmom*length(dist)/(areamm2));%% nAm/mm2
        if k==Npatchiter,
            fprintf('\n In last iteration...setting up patch %d with  %3.2f nAm , FWHM %3.2fmm, mean moment density %3.2f pAm/mm2, peak momemnt density %3.2f pAm/mm2 \n',j,nAm(j),smoothm(j)*1000,mompermm2(j)*1000,peakmompermm2(j)*1000);
        end;
        
    end; % for j
    
    fprintf('Prior %d. Average Mean (over dist of FWHM from centre) moment density %3.2f, sd %3.2f pAm/mm2\n',k,mean(mompermm2)*1000,std(mompermm2)*1000);
    fprintf('Prior %d. Average Peak (max vertex) moment density %3.2f, sd %3.2f pAm/mm2\n',k,mean(peakmompermm2)*1000,std(peakmompermm2)*1000);
    
    
    % NOW MAYBE WRITE A NEW PATCH FILE
    
    idnum=randi(1e6);
    priorfname=[priordir filesep sprintf('prior%d.mat',idnum)];
    fprintf('Saving %s\n',priorfname);
    F=[]; % no associated free energy value
    allpriornames=strvcat(allpriornames,priorfname);
    save(priorfname,'Qp','Qe','UL','F', spm_get_defaults('mat.format'));
end; % for k


%==========================================================================
function [q,dist,useind,areapervertex] = gauss_patch(M,i,FWHM,q)
%function [q,dist,useind,dist2fwhm,useind2]= gauss_patch(M,i,FWHM,q)


order=1;

sigma=FWHM./2.355;
sigma2=sigma^2;

d = spm_mesh_geodesic(M,i-1,order);

useind2=find(d<10/1000); %% a much larger area so that for small FWHM the per vertex area can be calculated
dist2fwhm=d(useind2);
areapervertex=pi*(max(dist2fwhm/2).^2)/length(useind2);


useind=find(d<=FWHM); %% maybe only extend just over 2 sigma in any one direction (i.e. cover 95% interval)

if FWHM==0
    useind=i;
    dist=0;
    q(useind)=1; %% impulse for 0 FWHM
else
    dist=d(useind);
    q(useind)=exp(-(dist.^2)/(2*sigma2));
end
