function [xXa,Sessa,Ka,Pa,nscana,rowa] = spm_bch_tsampl(xX,Sess,K,P,nscan,row)
% Used only in bch mode for the case of non regular temporal sampling
% @(#)spm_bch_tsampl.m	2.3 Stephanie Rouquette 99/09/15

global batch_mat;
global iA;

[xXa,Sessa,Ka,Pa,nscana,rowa] = deal(xX,Sess,K,P,nscan,row);

if ~isempty(batch_mat)
	xXa.X = [];
	Pa = '';
	compt1 = 0;
	compt2 = 0;
	for s = 1:length(Sess)
		sample = spm_input('batch',batch_mat,{'model',iA},'remain',s);
		nscana(s) = length(sample);
		xXa.X = [xXa.X;xX.X(compt1+sample,:)];
		Pa = [Pa;P(compt1+sample,:)];
		rowa{s} = (compt2+(1:nscana(s)))';
		Sessa{s}.rowa = rowa{s};
		Ka{s}.rowa = rowa{s};	
		compt1 = compt1+spm_input('batch',batch_mat,{'model',iA},'nscans',s);
		compt2 = compt2+nscan(s);
	end % for s = 1:length(Sess)
end


