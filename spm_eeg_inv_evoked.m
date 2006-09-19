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
% $Id: spm_eeg_inv_evoked.m 630 2006-09-19 15:21:04Z karl $


if D.events.Ntypes ~= D.Nevents
    error(sprintf('Evoked data required\n'));
end
try
    val = D.val;
catch
    val = length(D.inv);
end

qe = length(Qe);
qp = length(Qp);

% LOAD GAIN MATRIX
variabl = load(D.inv{val}.forward.gainmat);
name    = fieldnames(variabl);
G       = getfield(variabl, name{1});


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

woi      = D.inv{val}.inverse.woi;
contrast = D.inv{val}.inverse.contrast;
woi(1)   = round(woi(1)*(D.Radc/1000)) + D.events.start + 1;
woi(2)   = round(woi(2)*(D.Radc/1000)) + D.events.start + 1;
Nsens    = size(G,1);
Nsour    = size(G,2);

% get data
%--------------------------------------------------------------------------
y = sparse(0);
Y = sparse(0);
j = D.channels.eeg;
k = woi(1):woi(2);
t = D.events.start:D.events.stop;
for i = 1:D.Nevents
    l = find(D.events.code == D.events.types(i));
    y = y + contrast(i)*squeeze(D.data(j,k,l));
    Y = Y + contrast(i)*squeeze(D.data(j,t,l));
end

% DATA VARIANCE PARTITIONING (Call for ReML)
%==========================================================================
Yscal = 1/norm(Y,1);   % change the units of the data
Gscal = 1/norm(G,1);   % change the units of the sources

Y     = Yscal*Y;
y     = Yscal*y;
G     = Gscal*G;
YY    = Y*Y';
X     = ones(Nsens,1);
for i = 1:qe
    Qe{i}   = Qe{i}/norm(Qe{i},1);
end
for i = 1:qp
    Qp{i}   = Qp{i}/norm(Qp{i},1);
    GQpG{i} = G*Qp{i}*G';
end
Q      = {Qe{:} GQpG{:}};
Nbins  = size(Y,2);

% ReML (using all the data - Y)
%--------------------------------------------------------------------------
[C,h,Ph,F] = spm_reml(YY,X,Q,Nbins,0);

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


% COMPUTE THE SOURCE MAP OPERATOR
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

% compute MPA for time-windowed data (Y = y)
%--------------------------------------------------------------------------
Y     = y;
Nbins = size(Y,2);
J     = MAP*Y;

% Compute the proportion explained variance
%--------------------------------------------------------------------------
TotVar = sum(var(Y,  0,2));
ExpVar = sum(var(G*J,0,2));
R2     = 100*(ExpVar/TotVar)

% and rescale
%--------------------------------------------------------------------------
J      = Gscal*J/Yscal;

clear Y
if (D.inv{val}.mesh.Ctx_Nv - D.inv{val}.inverse.dim)
    Jfull       = sparse(D.inv{val}.mesh.Ctx_Nv,Nbins);
    Jfull(Is,:) = J;
    J           = Jfull;
    clear Jfull
end


% SAVE RESULTS
[pth,nam,ext] = spm_fileparts(D.fname);
woi           = D.inv{val}.inverse.woi;
D.inv{val}.inverse.resfile = [nam '_remlmat_' num2str(woi(1)) '_' num2str(woi(2)) 'ms_evoked.mat'];
Existence = exist(D.inv{val}.inverse.resfile);
if Existence == 2
    i = 2;
    while Existence == 2;
        D.inv{val}.inverse.resfile = [nam '_remlmat_' num2str(woi(1)) '_' num2str(woi(2)) 'ms_evoked_' num2str(i) '.mat'];
        Existence = exist(D.inv{val}.inverse.resfile);
        i = i + 1;
    end
end
D.inv{val}.inverse.LogEv   = F;

if spm_matlab_version_chk('7.1') >=0
    save(fullfile(pth,D.inv{val}.inverse.resfile), '-V6', 'C','h','Ph','F','J');
else
    save(fullfile(pth,D.inv{val}.inverse.resfile), 'C','h','Ph','F','J');
end

if spm_matlab_version_chk('7.1') >= 0
    save(fullfile(D.path,D.fname), '-V6','D');
else
    save(fullfile(D.path,D.fname), 'D');
end


% temporary visualisation
%--------------------------------------------------------------------------
[Finter,Fgraph] = spm('FnUIsetup','Stats: Results');
figure(Fgraph);
colormap('pink')

AvJev = abs(mean(J,2));
load(D.inv{D.val}.mesh.tess_ctx);

subplot(3,1,2)
patch('Vertices',vert,'Faces',face,'FaceVertexCData',AvJev,'FaceColor','flat');
view(-90,0);
shading interp
axis image
title('Averaged activity over the whole time window');

subplot(3,2,1)
patch('Vertices',vert,'Faces',face,'FaceVertexCData',AvJev,'FaceColor','flat');
view(-90,0);
shading interp
axis image
subplot(3,2,2)
patch('Vertices',vert,'Faces',face,'FaceVertexCData',AvJev,'FaceColor','flat');
view(90,0);
shading interp
axis image
subplot(3,2,5)
patch('Vertices',vert,'Faces',face,'FaceVertexCData',AvJev,'FaceColor','flat');
view(0,90);
shading interp
axis image
subplot(3,2,6)
patch('Vertices',vert,'Faces',face,'FaceVertexCData',AvJev,'FaceColor','flat');
view(0,-90);
shading interp
axis image

return



