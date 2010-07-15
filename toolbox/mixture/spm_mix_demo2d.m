
% Demonstrate use of spm_mix
m_true=3;
m_model=5;
d=2;
N=200;

disp('2-dimensional data;');
disp(sprintf('%d-cluster data and %d-cluster model',m_true,m_model));

% Use data from file
load yrep
figure
plot(y(:,1),y(:,2),'x');
hold on
   
vbmix=spm_mix(y,m_model);

for i=1:m_model,
   plot(vbmix.state(i).m(1),vbmix.state(i).m(2),'rx');
end
hold on
spm_mix_plot2d(vbmix,[-2 12 -2 12],1,'r',0.4,0.5);
set(gca,'FontSize',18);


