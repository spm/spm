function spm_eeg_inv_render(J,tess_ctx)

%=======================================================================
% Temproary display of source activity
%
% FORMAT spm_eeg_inv_render(J,tess_ctx)
% Input:
% J		    - conditional mean
% tess_ctx  - filename of 'vert' and 'face' file
%
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_render.m 630 2006-09-19 15:21:04Z karl $


% temporary visualisation
%--------------------------------------------------------------------------
J     = abs(mean(J,2));
load(tess_ctx);

patch('Vertices',vert,'Faces',face,'FaceVertexCData',J,'FaceColor','flat');
%view(0,0);
shading interp
axis image

return

% temporary visualisation
%--------------------------------------------------------------------------
[Finter,Fgraph] = spm('FnUIsetup','Stats: Results');
figure(Fgraph);
colormap('jet')

subplot(2,2,1)
patch('Vertices',vert,'Faces',face,'FaceVertexCData',J,'FaceColor','flat');
view(-90,0);
shading interp
axis image
subplot(2,2,2)
patch('Vertices',vert,'Faces',face,'FaceVertexCData',J,'FaceColor','flat');
view(90,0);
shading interp
axis image
subplot(2,2,3)
patch('Vertices',vert,'Faces',face,'FaceVertexCData',J,'FaceColor','flat');
view(0,90);
shading interp
axis image
subplot(2,2,4)
patch('Vertices',vert,'Faces',face,'FaceVertexCData',J,'FaceColor','flat');
view(0,-90);
shading interp
axis image

return



