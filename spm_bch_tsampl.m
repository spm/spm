function [xXa,Sessa,Ka,Pa,nscana,row] = spm_bch_tsampl(xX,Sess,K,P,nscan)

global batch_mat;
global iA;

[xXa,Sessa,Ka,Pa,nscana] = deal(xX,Sess,K,P,nscan);

if ~strcmp(batch_mat,'')
	xXa.X = [];
	Pa = '';
	compt1 = 0;
	compt2 = 0;
	for s = 1:length(Sess)
		sample = spm_input_b('batch',batch_mat,{'model',iA},'remain',s);
		nscana(s) = length(sample);
		xXa.X = [xXa.X;xX.X(compt1+sample,:)];
		Pa = [Pa;P(compt1+sample,:)];
		row{s} = (compt2+(1:nscana(s)))';
		Sessa{s}.row = row{s};
		Ka{s}.row = row{s};	
		compt1 = compt1+spm_input_b('batch',batch_mat,{'model',iA},'nscans',s);
		compt2 = compt2+nscan(s);
	end % for s = 1:length(Sess)
end


