function [xXa,Sessa,Ka,Pa,nscana,rowa] = spm_bch_tsampl(xX,Sess,K,P,nscan,row)
% Used only in bch mode for the case of non regular temporal sampling
% %W%	2.3 Stephanie Rouquette %E%

% Programmers Guide
% Batch system implemented on this routine. See spm_bch.man
% If inputs are modified in this routine, try to modify spm_bch.man
% and spm_bch_bchmat (if necessary) accordingly. 
% Calls to spm_input in this routine use the BCH gobal variable.  
%    BCH.bch_mat 
%    BCH.index0  = {'normalize',index_of_Analysis};
% or 
%    BCH.index0  = {'normalisation',index_of_Analysis}; (when
%                   spm_spn3d is launched for edit_defaults 
%
%_______________________________________________________________________

global BCH;

[xXa,Sessa,Ka,Pa,nscana,rowa] = deal(xX,Sess,K,P,nscan,row);

if ~isempty(BCH)
	xXa.X = [];
	Pa = '';
	compt1 = 0;
	compt2 = 0;
	for s = 1:length(Sess)
		sample = spm_input('batch',{},'remain',s);
		nscana(s) = length(sample);
		xXa.X = [xXa.X;xX.X(compt1+sample,:)];
		Pa = [Pa;P(compt1+sample,:)];
		rowa{s} = (compt2+(1:nscana(s)))';
		Sessa{s}.rowa = rowa{s};
		Ka{s}.rowa = rowa{s};	
		compt1 = compt1+spm_input('batch',{},'nscans',s);
		compt2 = compt2+nscan(s);
	end % for s = 1:length(Sess)
end


