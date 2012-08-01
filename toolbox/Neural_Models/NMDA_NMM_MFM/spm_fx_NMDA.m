function [f] = spm_fx_NMDA(x_V,x_G,P,M)

VN = 60;
%% Pyamidal Cells & interneuron NMDA receptos

mag_block = 1/(1 + 0.2*exp(-0.062*(exp(P.scale_NMDA)*x_V)));
f =  (x_G*(VN - x_V)*mag_block);
       