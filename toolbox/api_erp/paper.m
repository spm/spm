clear
clear functions

% empirical analyses; MMN
%--------------------------------------------------------------
load DCMP3                        % None
spm_dcm_erp(DCM)

load DCMP3F                       % Forward
spm_dcm_erp(DCM)

load DCMP3B                       % Backward
spm_dcm_erp(DCM)

load DCMP3FB                      % Forward and backward
spm_dcm_erp(DCM)

load DCMP3FBL                     % Forward, backward and Lateral
spm_dcm_erp(DCM)

% empirical analyses; Faces
%--------------------------------------------------------------
load DCMFace
spm_dcm_erp(DCM)


% simulations
%===========================================================================

% load target analysis
%---------------------------------------------------------------------------
load DCMP3FB
m     = 1;
n     = length(DCM.Qp.T);

% Change in B{1}
%---------------------------------------------------------------------------
B     = linspace(log(1/2),log(2),16);
EP    = {};
CP    = {};
Q     = spm_erp_pack(diag(DCM.Cp),m,n);
Q     = Q.B{1};
[j k] = find(Q == min(min(Q(find(Q)))));
xY    = DCM.xY;
e     = xY.y;
e(:)  = sqrt(DCM.Ce)*randn(length(DCM.Ce),1);
for i = 1:length(B)
    Ep           = spm_erp_pack(DCM.Ep,m,n);
    Ep.B{1}(j,k) = B(i);
    Ep           = spm_vec(Ep);
    y            = feval(DCM.M.IS,Ep,DCM.M,DCM.xU);
    xY.y         = y + e;
    [Ep,Cp]      = spm_nlsi_GN(DCM.M,DCM.xU,xY);
    EP{i}        = Ep;
    CP{i}        = Cp;
end

save paper_B

% Change in error
%---------------------------------------------------------------------------
B     = linspace(1/2,2,16);
EP    = {};
CP    = {};
y     = feval(DCM.M.IS,DCM.Ep,DCM.M,DCM.xU);
xY    = DCM.xY;
for i = 1:length(B)
    xY.y         = y + e*B(i);
    [Ep,Cp]      = spm_nlsi_GN(DCM.M,DCM.xU,xY);
    EP{i}        = Ep;
    CP{i}        = Cp;
end

save paper_E

% Graphics
%===========================================================================

% BMS - Graphics
%---------------------------------------------------------------------------
load DCMP3
F0 = DCM.F;
load DCMP3F
F(1) = DCM.F;
load DCMP3B
F(2) = DCM.F;
load DCMP3FB
F(3) = DCM.F;
load DCMP3FBL
F(4) = DCM.F;

F    = F - F0;                           % reference to 'no changes'
subplot(2,1,1)
bar(F)
axis square

% Change in B - Graphics
%---------------------------------------------------------------------------
load paper_B
l     = spm_invNcdf(.95);

% target parameter
%---------------------------------------------------------------------------
for i = 1:length(B)
    
    Qp   = spm_erp_pack(EP{i},m,n);
    E(i) = Qp.B{1}(j,k);
    Qp   = spm_erp_pack(diag(CP{i}),m,n);
    C(i) = l*sqrt(Qp.B{1}(j,k));
    
end
subplot(2,1,1)
fill([B fliplr(B)],[(E - C) fliplr(E + C)],...
           [1 1 1]*.9,'EdgeColor',[1 1 1])
hold on
plot(B,E,B,B,':')
hold off
axis square

% non-target parameter
%---------------------------------------------------------------------------
j     = 2;
k     = 4;
Q     = DCM.Qp.B{1}(j,k);
for i = 1:length(B)
    
    Qp   = spm_erp_pack(EP{i},m,n);
    E(i) = Qp.B{1}(j,k);
    Qp   = spm_erp_pack(diag(CP{i}),m,n);
    C(i) = l*sqrt(Qp.B{1}(j,k));
    
end
subplot(2,1,2)
fill([B fliplr(B)],[(E - C) fliplr(E + C)],...
         [1 1 1]*.9,'EdgeColor',[1 1 1])
hold on
plot(B,E,[B(1) B(end)],[Q Q],':')
hold off
axis square


% Change in E - Graphics
%---------------------------------------------------------------------------
load paper_E
l     = spm_invNcdf(.95);
Q     = DCM.Qp.B{1}(j,k);
for i = 1:length(B)
    
    Qp   = spm_erp_pack(EP{i},m,n);
    E(i) = Qp.B{1}(j,k);
    Qp   = spm_erp_pack(diag(CP{i}),m,n);
    C(i) = l*sqrt(Qp.B{1}(j,k));
    
end

subplot(2,1,1)
fill([B fliplr(B)],[(E - C) fliplr(E + C)],[1 1 1]*.9,'EdgeColor',[1 1 1])
hold on
plot(B,E,[B(1) B(end)],[Q Q],':')
hold off
axis square

