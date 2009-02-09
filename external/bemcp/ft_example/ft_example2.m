% create volume conductor starting from unit sphere
[pnt, tri] = icosahedron162;

vol = [];
vol.cond = [1 1/80 1];
vol.source = 1; % index of source compartment
vol.skin   = 3; % index of skin surface
% brain
vol.bnd(1).pnt = pnt*88;
vol.bnd(1).tri = tri;
% skull
vol.bnd(2).pnt = pnt*92;
vol.bnd(2).tri = tri;
% skin
vol.bnd(3).pnt = pnt*100;
vol.bnd(3).tri = tri;

% create the BEM system matrix
cfg = [];
cfg.method = 'bemcp';
vol1 = prepare_bemmodel(cfg, vol);


% create some random electrodes
pnt = randn(200,3);
pnt = pnt(pnt(:,3)>0, :);  % only those on the upper half
sens = [];
for i=1:size(pnt,1)
  sens.pnt(i,:) = pnt(i,:) / norm(pnt(i,:)); % scale towards the skin surface
  sens.label{i} = sprintf('%02d', i);
end

% prepare the sensor array and volume conduction, i.e. set up the linear
% interpolation from vertices to electrodes
[vol2, sens] = prepare_vol_sens(vol1, sens);

lf = compute_leadfield([0 0 50], sens, vol2);

figure; triplot(sens.pnt, [], lf(:,1));
figure; triplot(sens.pnt, [], lf(:,2));
figure; triplot(sens.pnt, [], lf(:,3));
