
% DCM for phase coupling: 
% Finger movement demo eg. Haken et al. 1985
%
% Sin-terms in PIF are (m=1) 1st or (m=2) 2nd order
%
% - if 1st order, sync is stable
% - if 2nd order, anti-sync (and sync) are stable

sim.m=1; % first order PIF
sim.do_plot=1;
sim.Nt=2; % Number of trials
sim.noise_dev=0.02;
sim.init='full'; % Uniform dist of phase diff between 0 and 2pi

[yy,P,M,U] = gen_finger (sim);

DCM = dcm_fit_finger (yy,M,U,sim.m);

disp(' ');
disp(' ');
for i=1:sim.m,
    disp(sprintf('Fourier order %d',i));
    disp('True parameters');
    disp(P.As(:,:,i));
    disp('DCM Estimates:');
    disp(abs(DCM.Ep.As(:,:,i)));
end

disp('GLM/EMA estimates:');
[A,fint] = glm_phi (yy,U.dt,M.fb);
disp(A);

Nr=length(M.x);
for j=1:Nr,
    figure
    for i=1:sim.Nt,
        subplot(sim.Nt,1,i);
        plot(U.tims,sin(yy{i}(:,j)));
        hold on
        plot(U.tims,sin(DCM.y{i}(:,j)),'r');
        legend('Data','DCM fit');
        if i==1
            title(sprintf('Region %d',j));
        end
    end
end




