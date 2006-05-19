function D = spm_eeg_inv_evoked(D,Qe,Qp)

%=======================================================================
% Inverse Solution for evoked activity
%
% FORMAT D = spm_eeg_inv_inverse_ui(D,Qe,Qp)
% Input:
% D		    - input data structure
% Qe        - noise covariance structure
% Qp        - prior source covaraince structure
% Output:
% D			- same data struct including the inverse solution files and variables
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_evoked.m 539 2006-05-19 17:59:30Z Darren $


if D.events.Ntypes ~= D.Nevents
    error(sprintf('Evoked data required\n'));
end
val = length(D.inv);

qe = length(Qe);
qp = length(Qp);


% LOAD GAIN MATRIX
variabl = load(D.inv{val}.forward.gainmat);
name    = fieldnames(variabl);
G       = getfield(variabl , name{1});


% SOURCE SPACE DIMENSION REDUCTION
if (D.inv{val}.mesh.Ctx_Nv - D.inv{val}.inverse.dim)
    Msize = D.inv{val}.inverse.dim
    if isempty(D.inv{val}.inverse.priors.level2{3}.filename)
        error(sprintf('Multivariate Source Prelocalisation has not been performed\n'));
    end
    variabl   = load(D.inv{val}.inverse.priors.level2{3}.filename);
    name      = fieldnames(variabl);
    APM       = getfield(variabl , 'APM');
    [APMs,Is] = sort(-APM);
    APMs      = -APMs;
    Is        = Is(1:Msize);
    G         = G(:,Is);
    GQpG      = {};
    for i = 1:qp
        Qp{i} = Qp{i}(Is,Is);
    end
end

woi        = D.inv{val}.inverse.woi;
contrast   = D.inv{val}.inverse.contrast;
Nsens      = D.Nchannels;
Nsour      = size(G,2);

woi(1) = round(woi(1)*(D.Radc/1000)) + D.events.start + 1;
woi(2) = round(woi(2)*(D.Radc/1000)) + D.events.start + 1;

Y     = sparse(zeros(size(D.data,1),(woi(2)-woi(1)+1)));
for i = 1:D.Nevents
    k = find(D.events.code == D.events.types(i));
    Y = Y + contrast(i)*squeeze(D.data(:,woi(1):woi(2),k));
end

% DATA VARIANCE PARTITIONING (Call for ReML)
ExpScal     = max(Y(:));
Scal        = floor( log10(ExpScal) );
Scal        = 10^(-Scal);

Y     = Scal*Y;
G     = Scal*G;
YY    = Y*Y';
X     = ones(Nsens,1);
for i = 1:qp
    GQpG{i} = G*Qp{i}*G';
end
Q      = {Qe{:} GQpG{:}};
Nsamp  = woi(2) - woi(1) + 1;

[C,h,Ph,F] = spm_reml(YY,X,Q,Nsamp,1);
clear YY Q


% DISPLAY ReML OUTCOME
disp(['Model Log-Evidence:  ' num2str(F)]);
DisplayString = ['Noise covariance hyperparameters:  '];
for i = 1:qe
    DisplayString = [DisplayString '  ' num2str(h(i))];
end
disp(DisplayString);
DisplayString = ['Source covariance hyperparameters:  '];
for i = 1:qp
    DisplayString = [DisplayString '  ' num2str(h(qe + i))];
end
disp(DisplayString);


% COMPUTE THE SOURCE MAP ESTIMATE
Ce = sparse(Nsens,Nsens);
for i = 1:qe
   Ce = Ce + h(i)*Qe{i};
end
clear Qe GQpG
Cp = sparse(Nsour,Nsour);
for i = 1:qp
    Cp = Cp + h(qe + i)*Qp{i};
end
clear Qp

CpG  = Cp*G';
GCpG = G*CpG;
MAP  = CpG*inv(GCpG + Ce);
clear CpG GCpG Ce Cp

J = MAP*Y;


% Compute the explained variance
TotVar = norm(Y - mean(Y,1)*ones(Nsamp,1),'fro')^2;
ExpVar = TotVar - norm(Y - G*J,'fro')^2;
R2 = 100*(ExpVar/TotVar)


clear Y
if (D.inv{val}.mesh.Ctx_Nv - D.inv{val}.inverse.dim)
    Jfull = sparse(D.inv{val}.mesh.Ctx_Nv,Nsamp);
    Jfull(Is,:) = J;
    J = Jfull;
    clear Jfull
end


% SAVE RESULTS
[pth,nam,ext] = spm_fileparts(D.fname);
woi           = D.inv{val}.inverse.woi;
Ntime         = clock;
D.inv{val}.inverse.resfile = [nam '_remlmat_' num2str(woi(1)) '_' num2str(woi(2)) 'ms_evoked_' num2str(Ntime(4)) 'H' num2str(Ntime(5)) '.mat'];
D.inv{val}.inverse.LogEv   = F;

if spm_matlab_version_chk('7') >=0
    save(fullfile(pth,D.inv{val}.inverse.resfile), '-V6', 'C','h','Ph','F','J');
else
    save(fullfile(pth,D.inv{val}.inverse.resfile), 'C','h','Ph','F','J');
end

if spm_matlab_version_chk('7') >= 0
    save(fullfile(D.path,D.fname), '-V6','D');
else
    save(fullfile(D.path,D.fname), 'D');
end

% Temporary Visualization
AvJev = mean(J')';
load(D.inv{val}.mesh.tess_ctx);
colormap jet
axis off
patch('Vertices',vert,'Faces',face,'FaceVertexCData',AvJev,'FaceColor','flat');
view(-90,0);
shading interp
colorbar
title('Averaged activity over the whole time window');
cameramenu
