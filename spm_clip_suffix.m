function C = spm_clip_suffix(P)
% clips suffix from filename
% FORMAT T = spm_clip_suffix(P)
% P - filename
%____________________________________________________________________________
%
% spm_clip_suffix returns the file name P with the any recognized
% suffix removed.  Any spaces are removed from the filename.
% 
% Recognized sufficies are ".img" ".hdr" ".air"
%
% from spm_clip_suffix.m, v1.2, UPMC PET Modified SPM (6/24/95)
% %W% %E%

% Recognized sufficies... all must be same length
Suf = [	'.img';
	'.hdr';
	'.air'	];

% Remove spaces
P     = deblank(P);

% See if suffix of P matches Suf
q     = length(P);
Pt    = P(q-size(Suf,2)+1:q);
if find(sum((Pt(ones(1,size(Suf,1)),:) ~= Suf)') == 0) ~= []
	C = P(1:(q - 4));
else
	C = P;
end
