function [a,b,c,d] = tf2ss(num, den)
%TF2SS  Transfer function to state-space conversion.
%   [A,B,C,D] = TF2SS(NUM,DEN)  calculates the state-space 
%   representation:
%       .
%       x = Ax + Bu
%       y = Cx + Du
%
%   of the system:
%               NUM(s) 
%       H(s) = --------
%               DEN(s)
%
%   from a single input.  Vector DEN must contain the coefficients of
%   the denominator in descending powers of s.  Matrix NUM must 
%   contain the numerator coefficients with as many rows as there are
%   outputs y.  The A,B,C,D matrices are returned in controller 
%   canonical form.  This calculation also works for discrete systems.
%
%   For discrete-time transfer functions, it is highly recommended to
%   make the length of the numerator and denominator equal to ensure 
%   correct results.  You can do this using the function EQTFLENGTH in
%   the Signal Processing Toolbox.  However, this function only handles
%   single-input single-output systems.
%
%   See also TF2ZP, SS2TF, ZP2SS, ZP2TF.

%   J.N. Little 3-24-85
%   Copyright 1984-2003 The MathWorks, Inc.
%   $Revision: 1.22.4.1 $  $Date: 2003/05/01 20:42:48 $
%   Latest revision 4-29-89 JNL, 7-29-96 PG

[mnum,nnum] = size(num);
[mden,n] = size(den);
% Check for null systems
if  (n == 0 && nnum == 0), a=[]; b=[]; c=[]; d=[]; return, end

if min(mden,n)>1,
   % Error out if DEN is an array
   error('MATLAB:tf2ss:NeedRowDenom', 'Denominator must be a row vector.');
elseif mden>1,
   % Transpose DEN when a column vector
   den = den.';
end

% Strip leading zeros from denominator
inz = find(den ~= 0);
den = den(inz(1):end);
[mden,n] = size(den);

% Check for proper numerator
if nnum > n
    % Try to strip leading zeros to make proper
    if (all(all(num(:,1:(nnum-n)) == 0)))
        num = num(:,(nnum-n+1):nnum);
        [mnum,nnum] = size(num);
    else
        error('MATLAB:tf2ss:DenomInvalidOrder',...
              ['Order of denominator must be greater than or equal to ',...
          'order of numerator.']);
    end
end

% Pad numerator with leading zeros, to make it have the same number of
% columns as the denominator, and normalize it to den(1)
num = [zeros(mnum,n-nnum) num]./den(1);

% Do the D-matrix first
if length(num)
    d = num(:,1);
else
    d = [];
end

% Handle special constant case:
if n == 1
    a = [];
    b = [];
    c = [];
    return
end

% Now do the rest, starting by normalizing den to den(1),
den = den(2:n) ./ den(1);
a = [-den; eye(n-2,n-1)];
b = eye(n-1,1);
if mnum > 0
    c = num(:,2:n) - num(:,1) * den;
else
    c = [];
end

