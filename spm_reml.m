function [Ce,h] = spm_reml(Cy,X,Q);
% REML estimation of covariance components from Cov{y}
% FORMAT [Ce,h] = spm_reml(Cy,X,Q);
%
% Cy  - (m x m) data covariance matrix y*y'
% X   - (m x p) design matrix
% Q   - {1 x q} covariance constraints
%
% Ce  - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) hyperparameters
%___________________________________________________________________________
% %W% John Ashburner, Karl Friston %E%

q = length(Q);
o = size(Q{1},1);

% Find any block structure in the data, and divide it into cells of QQ.
%---------------------------------------------------------------------------
msk = logical(zeros(size(Q{1},1),length(Q)));
for i=1:length(Q), msk(:,i)=any(Q{i})'; end;
lkp = 1:length(Q);
QQ  = {};
ii  = 1;
while length(Q),
	msk0 = msk(:,1);
	for j=1:size(msk,2),
		for i=1:size(msk,2),
			if any(msk(:,i)&msk0),
				msk0 = msk0|msk(:,i);
			end;
		end;
	end;
	sel = logical(zeros(size(Q,2),1));
	for i=1:size(msk,2),
		sel(i) = any(msk(:,i)&msk0);
	end;
	if any(sel),
		QQ{ii}.Q   = Q(sel);
		QQ{ii}.lkp = lkp(sel);
		QQ{ii}.msk = msk0;
		for i=1:length(QQ{ii}.Q),
			QQ{ii}.Q{i} = QQ{ii}.Q{i}(msk0,msk0);
		end;
		Q   = Q(~sel);
		lkp = lkp(~sel);
		msk = msk(:,~sel);
		ii  = ii+1;
	end;
end;
clear Q;

% ensure X is not rank deficient
%---------------------------------------------------------------------------
X     = spm_svd(X);

% REML	objective function 2F = -r'*iCe*r - log(det(Ce)) - log(det(XiCeX));
%---------------------------------------------------------------------------
J     = zeros(q,q);
g     = zeros(q,1);
h     = zeros(q,1);
for ii=1:length(QQ),
	for i = 1:length(QQ{ii}.Q),
		h(QQ{ii}.lkp(i),1) = full(any(diag(QQ{ii}.Q{i})));
	end;
end;

for k = 1:32

	Cby = zeros(size(X,2));
	for ii=1:length(QQ),
		QQ{ii}.Ce = sparse(sum(QQ{ii}.msk),sum(QQ{ii}.msk));
		for i = 1:length(QQ{ii}.Q),
			QQ{ii}.Ce = QQ{ii}.Ce + h(QQ{ii}.lkp(i))*QQ{ii}.Q{i};
		end;
		QQ{ii}.iCe   = inv(full(QQ{ii}.Ce));
		QQ{ii}.iCeX  = QQ{ii}.iCe*X(QQ{ii}.msk,:);
		Cby          = Cby + X(QQ{ii}.msk,:)'*QQ{ii}.iCeX;
	end;
	Cby = inv(Cby);

	% 1st derivatives g = 2dF/dhi
	% 2nd derivatives J = 2ddF/didhj
	%------------------------------------------------------------------
	for ii=1:length(QQ),
		QQ{ii}.R = QQ{ii}.iCe - QQ{ii}.iCeX*Cby*QQ{ii}.iCeX';
		RRCyR    = QQ{ii}.R;
		R        = cell(length(QQ),1);
		for jj=1:length(QQ),
			if jj~=ii, R{jj} = -QQ{ii}.iCeX*Cby*QQ{jj}.iCeX';
			else,      R{jj} =  QQ{ii}.R; end;
		end;
		for jj=1:length(QQ),
			for kk=jj:length(QQ),
				tmp = R{jj}*Cy(QQ{jj}.msk,QQ{kk}.msk)*R{kk}';
				if jj==kk, RRCyR = RRCyR - tmp;
				else,      RRCyR = RRCyR - tmp - tmp'; end;
			end;
		end;

		for i = 1:length(QQ{ii}.Q),
			lkp       = QQ{ii}.lkp;
			g(lkp(i)) = sum(sum(RRCyR.*QQ{ii}.Q{i}));
			QC        = QQ{ii}.Q{i}*QQ{ii}.iCe;
			J(lkp(i),lkp(i)) = sum(sum(QC.*QC));
			for j  = (i+1):length(QQ{ii}.Q),
				J(lkp(i),lkp(j)) = sum(sum(QC.*(QQ{ii}.Q{j}*QQ{ii}.iCe)));
				J(lkp(j),lkp(i)) = J(lkp(i),lkp(j));
			end;
		end;
	end

	% Newton like step
	%------------------------------------------------------------------
	dh    = -J\g(:);
	h     = h + dh;

	% Convergence
	%===================================================================
	w     = dh'*dh;
	fprintf('%-30s: %i %30s%e\n','  ReML Iteration',k,'...',full(w));
	if w < 1e-6, break, end
end

% Reconstitute Covariance Matrix
% --------------------------------------------------------------------------
Ce = sparse(o,o);
for ii=1:length(QQ),
	for jj=1:length(QQ{ii}.Q),
		Ce(QQ{ii}.msk,QQ{ii}.msk) = Ce(QQ{ii}.msk,QQ{ii}.msk) + h(QQ{ii}.lkp(jj))*QQ{ii}.Q{jj};
	end;
end;
