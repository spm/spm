function [Lc2z,Lc2t,Lt2z] = spm_lambda(df)
% Coefficient relating var(dR/dx), var(dt/dx) and var(dZ/dx)
% FORMAT [Lc2z,Lc2t,Lt2z] = spm_lambda(df)
% df   - degrees of freedom of the t field
%
% Lc2z - coefficient relating component to Gaussianized t fields
% Lc2t - coefficient relating component to t fields
% Lt2z - coefficient relating t fields  to Gaussianized t fields
%_______________________________________________________________________
%
%
% returns the coefficients relating var(dR/dx), var(dt/dx) and
% var(dZ/dx), where R are the component fields of a t-field t, of degrees
% of freedom, df and Z is the Gaissianized t field.  dZ/dx are the first
% spatial partial derivatives of Z.
%
% ref:  Worsley KJ et al (1992) JCBFM 12 p917
%       Holmes AP (1994) PhD Thesis p224 App.G: Smoothness of t-fields
%_______________________________________________________________________
% %W% Karl Friston %E%


% for high df
%-----------------------------------------------------------------------
if df > 32; df = round(df); end

% pdfs and Probability integral transform (PIT)
%-----------------------------------------------------------------------
dt    = 0.01;
t     = -12:dt:12;
pdf   = spm_Tpdf(t,df);
T     = spm_t2z(t,df);

% var(dt/dx) following transformation = var(dT(t)/dx)/var(dt/dx)
%-----------------------------------------------------------------------
Lc2z  = sum((t.^2 + df).^2/(df*(df - 1)).*pdf.^3./(spm_Npdf(T).^2)*dt);
Lc2t  = sum((t.^2 + df).^2/(df*(df - 1)).*pdf*dt);
Lt2z  = Lc2z/Lc2t;

% for low df
%-----------------------------------------------------------------------
LC2T  = [Inf Inf 20.98 4.37 2.66 2.09 1.80];
if df < length(LC2T); Lc2t = LC2T(df); end
