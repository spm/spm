function spm_eeg_invert_display(D,PST,Ndip)
% Displays conditional expectation of response (J)
% FORMAT spm_eeg_invert_display(D,PST,Ndip)
% FORMAT spm_eeg_invert_display(D,XYZ,Ndip)
% D    - 3D structure (ReML estimation of response (J) )
% PST  - peristimulus time (ms) - defaults to the PST of max abs(J)
%      - [Start Stop] (ms)     - invokes a movie of CSD
% XYZ  - dipole location of interest
%
% Ndip - number of dipole to display (default 512)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_eeg_invert_display.m 5219 2013-01-29 17:07:07Z spm $
 
% Number of dipoles to display
%==========================================================================
try, PST;  catch, PST  = [];  end
try, Ndip; catch, Ndip = 512; end
 
% get condition
%--------------------------------------------------------------------------
if ~isfield(D, 'con') || D.con == 0
    con = 1;
else
    con = D.con;
end
 
% D - SPM data structure
%==========================================================================
model = D.inv{D.val};

con   = min(con,length(model.inverse.J));

try
    disp(model.inverse);
catch
    warndlg('please invert model')
    return
end
 
% get solution and spatiotemporal basis
%--------------------------------------------------------------------------
J      = model.inverse.J;
T      = model.inverse.T;
Is     = model.inverse.Is;
pst    = model.inverse.pst;
R2     = model.inverse.R2;
F      = model.inverse.F;
Nd     = model.inverse.Nd;

Ndip   = min(Ndip,length(Is));
try
    VE = model.inverse.VE;
catch
    VE = R2;
end

 
% - project J onto pst
%--------------------------------------------------------------------------
J      = J{con}*T';
 
% display
%==========================================================================
Fgraph = spm_figure('GetWin','Graphics');
figure(Fgraph);clf
vert   = model.mesh.tess_mni.vert;
 
% movie
%--------------------------------------------------------------------------
if length(PST) == 2
    
    
    % get significant voxels
    %----------------------------------------------------------------------
    Nb     = 170;
    [i,j1] = min(abs(pst - PST(1)));
    [i,j2] = min(abs(pst - PST(2)));
    jt     = fix(linspace(j1,j2,Nb));
    J      = abs(J(:,jt));
    Z      = max(J,[],2);
    [T,j]  = sort(-Z);
    js     = j(1:Ndip);
    
    J      = J(js,:);
    XYZ    = vert(Is(js),:)';
    Jmax   = max(max(J));
    Jmin   = min(min(J));
    XYZ    = [XYZ [ones(1,Nb)*(-128); ([1:Nb] - 100); ones(1,Nb)*(-64)]];
    
    % MIP
    %----------------------------------------------------------------------
    subplot(2,1,1)
    for j = 1:Nb
        SCL      = ones(Nb,1)*Jmin;
        SCL(1:j) = Jmax;
        figure(Fgraph)
        spm_mip([J(:,j); SCL],XYZ,6);
        axis image
        title({sprintf('Response at %i most active voxels',Ndip),...
               sprintf('at %i ms (from %i to %i ms)',fix(pst(jt(j))),fix(pst(j1)),fix(pst(j2)))})
        drawnow
    end
    return
end
 
% maximum response at XYZ
%--------------------------------------------------------------------------
if length(PST) == 3
    [i,js] = min(sum([vert(Is,1) - PST(1), ...
                      vert(Is,2) - PST(2), ...
                      vert(Is,3) - PST(3)].^2,2));
    [i,jt] = max(abs(J(js,:)));
    
% maximum response at PST
%--------------------------------------------------------------------------
elseif length(PST) == 1
    [i,jt] = min(abs(pst - PST));
    [i,js] = max(abs(J(:,jt)));
    
% maximum response in space and time
%--------------------------------------------------------------------------
else
    [i,js] = max(max(abs(J),[],2));
    [i,jt] = max(abs(J(js,:)));
end
Jt    = J(js,:);                     % over time
Js    = J(:,jt);                     % over sources
PST   = round(pst(jt));
XYZ   = round(vert(Is(js),:));
Jmax  = abs(sparse(Is,1,Js,Nd,1));
 
% plot responses over time
%==========================================================================
subplot(2,1,1)
for i = 1:length(model.inverse.J)
    
    % i-th trial
    %----------------------------------------------------------------------
    J     = model.inverse.J{i};
    Jt    = J(js,:)*T';
    maxJ  = max(max(abs(Jt)));
    if i == con
        Color = [1 0 0];
    else
        Color = [1 1 1]*(1 - 1/4);
    end
    
    % plot confidence intervals is possible
    %----------------------------------------------------------------------
    try
        qC  = model.inverse.qC(js).*diag(model.inverse.qV)';
        ci  = 1.64*sqrt(abs(qC));
        plot(pst,Jt,pst,Jt + ci,':',pst,Jt - ci,':',...
            [PST PST],[-1 1]*maxJ,':',...
            'Color',Color)
        hold on
    catch
        plot(pst,Jt,[PST PST],[-1 1]*maxJ,':','Color',Color)
        hold on
    end
end
title({sprintf('estimated response - condition %d',con), ...
       sprintf('at %i, %i, %i mm',XYZ(1),XYZ(2),XYZ(3))})
xlabel('time  ms')
axis square
hold off
 
 
% Mean square responses over space
%==========================================================================
subplot(2,1,2)
[T,i]  = sort(-abs(Js));
i      = i(1:Ndip);
spm_mip(Jmax(Is(i)),vert(Is(i),:)',6);
axis image
 
% title
%--------------------------------------------------------------------------
try
    qC = model.inverse.qC(i).*model.inverse.qV(jt,jt);
    Z  = min(Jmax(Is(i))./sqrt(qC));
    PP = fix(100*(spm_Ncdf(Z)));
    title({sprintf('PPM at %i ms (%i percent confidence)',PST,PP), ...
           sprintf('%i dipoles',length(i)), ...
           sprintf('Percent variance explained %.2f (%.2f)',full(R2),full(VE)), ...
           sprintf('log-evidence = %.1f',full(F))})
catch
    title({sprintf('Responses at %i dipoles',length(i)), ...
           sprintf('Variance explained %.2f (percent)',full(R2)), ...
           sprintf('log-evidence = %.1f',full(F))})
end
drawnow
