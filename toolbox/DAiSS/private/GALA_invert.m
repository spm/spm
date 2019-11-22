function [J, S] = GALA_invert(BF,tS)

YTY   = sparse(0);                             

Nl = tS.Nl;

for p=1:Nl
    
    BFp = load(BF{p});
    
    % load data
    D{p} = spm_eeg_load(BFp.data.D);

    % load gain matrix
    modality='MEG'; % temporary - should be replaced by tS.modality
    L1 = BFp.sources.L.(modality);
    Nd = length(L1);
    Lp{p} = [];
    for i = 1:Nd
        lf = L1{i};     
        Lp{p} = cat(2, Lp{p}, lf);
    end
    % eliminate low SNR spatial modes
    U = spm_svd((Lp{p}*Lp{p}'),exp(-16));
    ULp{p} = U'*Lp{p};
    Nm(p) = size(ULp{p},1);
    Up{p} = U;

    trial = D{p}.condlist;      % get conditions list
    Nt = 2;                     % temporary invert only first one
    Nb = size(D{p},2);          % Number of time bins
    
    % get (spatially aligned) data
    YY = 0; N=0;
    for j = 1:Nt                                % loop over conditions
        Y{p,j}  = Up{p}'*D{p}(:,:,j);
        YY    = YY + Y{p,j}'*Y{p,j};
        N     = N + Nb;
    end

    YTY        = YTY + YY;
    
end

% Initialize matrices comtrolling covariance reparametrization
Is = 1:Nd;                          
for p=1:Nl
    Isp{p} = 1:Nd;                          
end
lIs = 1:Nd*Nl;                      
Vs = speye(Nd*Nl);                  

mNm = mean(Nm);
scale = 1/sqrt(trace(YTY)/(N*mNm));
YTY = YTY*(scale^2);

% temporal projector 
[V E]  = spm_svd(YTY,exp(-8));              % get temporal modes
E      = diag(E)/trace(YTY);                % normalise variance
Nr  = length(E);                            % number of temporal modes
S   = V(:,1:Nr);                            % temporal modes
VE  = sum(E(1:Nr));                         % variance explained
    
fprintf('Using %i temporal modes, ',Nr)
fprintf('accounting for %0.2f percent average variance\n',full(100*VE))

% scale and align data (spatial and temporary modes)
for p=1:Nl
    MY=[];
    for j = 1:Nt                % loop over Nt conditions
        MY{j}   = Y{p,j}*S*scale;
    end
    UYp{p} = spm_cat(MY);
end

% prepare long data and leadfields and channel noise components
UL=[]; Qe=[];
for p=1:Nl
    UL = blkdiag(UL,ULp{p});
    Qep{p} = Up{p}'*Up{p};
    Qe = full(blkdiag(Qe,Qep{p}));
end
UY = spm_cat(UYp');

Scale = sqrt(trace(UL*UL')/mNm);
UL = UL/Scale;

% prepare smoothing kernel - it could be done as in spm but this scheme
% applicable for FreeSurfer meshes
% -----------------------------------------------------------------------
radius = 4;

vert  = D{1}.inv{D{1}.val}.mesh.tess_mni.vert;
face  = D{1}.inv{D{1}.val}.mesh.tess_mni.face;

AA=[];
AAA     = spm_mesh_distmtx(struct('vertices',vert,'faces',face),0);
% find subsets of certain distance
AA{1}=speye(length(vert));
AA{2}=AAA;
B = AA{1}+AA{2};
A = B;
for i=2:radius-1
    in = find(AAA^i);
    Bn=sparse(length(vert),length(vert));
    Bn(in)=1;
    AA{i+1}=full(Bn-B);
    B=Bn;
end

y=[1 0.8 0.37 0.2];

QG = sparse(length(vert),length(vert));
for i=1:radius
    QG=QG + AA{i}*y(i);
end
QG =sparse(QG(:,:));
K    = sparse(QG.*(QG > exp(-8)));
clear AA AAA B Bn QG
% -----------------------------------------------------------------------

% covariance of the data
YY = UY*UY';

% noise covariance component
Qn{1} = Qe; % it's the first one for channel noise

% template for source-noise covariance
sQt = kron(speye(Nl),speye(Nd));

thresh = 0.5;       % threshold 

for it=1:tS.iter

    % starting from 2nd iteration add source-noise covariance matrices
    if it>1
        sQn = Vn*sQt*Vn; %Vn - noise-vertices matrix, defined below 
        Qn{it} = UL*sQn*UL'; 
    end

    % first ROI covariance matrix - strong correlation between subjects
    psQ1 = kron(ones(Nl),K);
%     psQ1 = kron(ones(Nl),speye(Nd));
    sQ1 = Vs*psQ1*Vs;
    Q1 = full(UL*sQ1*UL');
    
    % second ROI covariance matrix - subjects specific activity
    psQ2 = kron(speye(Nl),K);
%     psQ2 = kron(speye(Nl),speye(Nd));
    sQ2 = Vs*psQ2*Vs;
    Q2 = full(UL*sQ2*UL');
    
    Qs = {Q1 Q2};
%     Qs = {Q2};

    % ML estimation of hyperparameters
    [Cy,h,Ph,F] = spm_reml_sc(YY,[],[Qn Qs],1,-4,16);
    % take ROI components and discard noise components
    DD = h(end-1)*sQ1 + h(end)*sQ2;
%     DD = h(end)*sQ2;
    Cy = full(Cy);    

    M = DD*UL'/Cy;

    % long (all subjects) source activity
    J = M*UY;

    % goodness of fit
    SSR  = sum(var((UY - UL*J),0,2));
    SST  = sum(var( UY,0,2));
    R2   = 100*(SST - SSR)/SST;
    fprintf('\nPercent variance explained GALA %.2f\n',full(R2));
    fprintf('F = %.10f\n',full(F));
%     fprintf('L = %.10f\n',full(-spm_logdet(Cy)- trace(UY'*inv(Cy)*UY)));
%     fprintf('ML fit = %.10f\n',-trace(UY'*inv(Cy)*UY));
%     fprintf('ML det = %.10f\n\n',full(-spm_logdet(Cy))); 

    % idea - get vertices with max covariance both within and between
    % subjects. As the result they should not be the same (indices) for
    % different subjects, so correlated between subjects part of covarince
    % prior is overlapping but not identical patches

    % pay attention - it's to exclude disconected vertices from Jcov
    % it seems ssQ1 = spones(sQ1) works better than simple sQ1
    % may be because in Jcov.*sQ1 there is double attenuation of tails
    
    
    if it<tS.iter
        
        ssQ1 = spones(sQ1);

        % real calculation for big matrices
        J = full(J);

        lJcov = zeros(1,Nl*Nd);
        for i=1:Nd*Nl
            if any(J(i,:))
                Jcovi = J*J(i,:)';
                sJcovi = Jcovi.*ssQ1(:,i);
                lJcov(i) = squeeze(sum(sJcovi));
            end
        end

        % reparametrization of covariance matrices
        lIso = lIs;
        lIs = find(lJcov>quantile(lJcov,thresh));   % ROI vertices update
        lnIs = setdiff(lIso,lIs);                   % noise vertices update

        Vs = zeros(Nd*Nl,1);    
        Vn = zeros(Nd*Nl,1);    
        Vs(lIs) = 1;
        Vn(lnIs) = 1;
        Vs = spdiags(Vs,0,Nd*Nl,Nd*Nl);     % ROI vertices matrix update
        Vn = spdiags(Vn,0,Nd*Nl,Nd*Nl);     % noise vertices matrix update

        thresh = thresh+(1-thresh)/2;               % threshold update
    
    end
    
    % display
    
%     MM.vertices = vert;
%     MM.faces = face;
% 
%     ahp=[];
%     figure;
%     for p=1:Nl
%         lJcovp{p} = lJcov(1+(p-1)*Nd:Nd+(p-1)*Nd);
%         
%         ahp(p)=subplot(2,2,p);
%         srcs_disp1 = lJcovp{p}/max(abs(lJcovp{p}));
% 
%         cla; axis off
% 
%         fig1 = patch('vertices',MM.vertices,'faces',MM.faces,'FaceVertexCData',srcs_disp1');
% 
%         set(fig1,'FaceColor',[.5 .5 .5],'EdgeColor','none');
%         shading interp
%         lighting gouraud
%         %camlight
%         zoom off
%         lightangle(0,270);lightangle(270,0),lightangle(90,0),lightangle(0,45),lightangle(0,135);
%         material([.1 .1 .4 .5 .4]);
%         %view(140,15);
%         caxis([0 1])
%         colormap(jet);
% 
%     end
%     
%     hlink = linkprop(ahp, {'CameraPosition','CameraUpVector'});
%     key = 'graphics_linkprop';
%     setappdata(ahp(1),key,hlink); 
    
end

end



