function spm_eeg_invert_display(D)
% Diplays condiotnal expecation of response (J)
% FORMAT spm_eeg_invert_display(D)
% ReML estimation of response (J) 
%__________________________________________________________________________


% D - SPM data structure
%==========================================================================
model = D.inv{D.val};
try
    disp(model.inverse);
catch
    warndlg('please invert model')
    return
end

%--------------------------------------------------------------------------
J     = model.inverse.J;
qC    = model.inverse.qC;
Is    = model.inverse.Is;
Nd    = model.inverse.Nd;
Nt    = model.inverse.Nt;
pst   = model.inverse.pst;
R2    = model.inverse.R2;
F     = model.inverse.F;

% display
%==========================================================================
Fgraph   = spm_figure('GetWin','Graphics');
clf(Fgraph)
figure(Fgraph)
tess_ctx = model.mesh.tess_ctx;
load(model.mesh.tess_ctx,'vert')

% maximum response and confidence intervals
%--------------------------------------------------------------------------
[i j] = max(max(abs(J),[],2));
Jt    = J(j,:);                      % over time
ci    = Nt*sqrt(qC(j))*1.64;
[i j] = max(max(abs(J),[],1));
Js    = J(:,j);                      % over sources
Jmax  = abs(sparse(Is,1,Js,Nd,1));

% maximum response - space
%--------------------------------------------------------------------------
subplot(2,3,1)
spm_eeg_inv_render(Jmax,tess_ctx)
view([180 -90])
title('estimated response')

% maximum response - time
%--------------------------------------------------------------------------
subplot(2,2,2)
plot(pst,Jt,pst,Jt + ci,':',pst,Jt - ci,':')
xlabel('time  ms')
title('estimated response (95% intervals)')
axis square

% PPM
%==========================================================================
subplot(2,1,2)
Z     = abs(Js)./(Nt*sqrt(qC));
i     = find(Z > 1.64);

spm_mip(Jmax(Is(i)),vert(Is(i),:)',6);
axis image
title({'PPM estimated response (95% confidence)'...
       sprintf('Variance explained %.2f (percent)',R2)...
       sprintf('log-evidence = %.3f',F)})
drawnow

