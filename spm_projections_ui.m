function spm_projections_ui
% used to review results of statistical analysis (SPM{Z})
% FORMAT spm_projections_ui
%___________________________________________________________________________
%
% spm_projections_ui allows the SPM{Z} created by spm_spm.m to be re-displayed
% and characterized in terms of regionally significant effects.  Note that
% the threshold does not have to be the same as in the original analysis
%
% Multiple [orthogonal] contrasts can be specified to produce a SPM{Z} that
% reflects the significance of two or more effects:-
%
% specifying a vector for the contrasts causes the second (and ensuing)
% contrasts to mask the first.  Non-orthogonal compounds are allowed.
% The ensuing voxels reach criteria (at the uncorrected threshold
% specified) for all the contrasts.  The statistic constituting the
% SPM{Z} is the mean of the n component Z values divided by sqrt(n).
% Given the orthogonality constraint above, this statistic is itself
% Gaussian.  This use of spm_projections_ui.m is useful for testing
% multiple hypothese simultaneously.  For example the conjunction
% of activations in several task-pairs or (in multifactorial designs)
% the analysis of interaction effects in (and only in) areas subject to a
% main effect.  The resulting p values are generally so small that one
% can forgo a correction for multiple comparisons.
%
% see spm_projections.m for further details
%
%__________________________________________________________________________
% %W% %E%


% get CWD and results
%---------------------------------------------------------------------------
global CWD
set(2,'Name','SPM{Z} projections')

tmp   = spm_get(1,'.mat','select SPMt.mat for analysis','SPMt');
CWD   = strrep(tmp,'/SPMt.mat','');
K     = [];

load([CWD,'/SPM'])
load([CWD,'/XYZ'])
load([CWD,'/SPMt'])

% get contrast[s]
%---------------------------------------------------------------------------
i    = 0;
while any(i < 1 | i > size(CONTRAST,1))
	i = spm_input(sprintf('contrast[s] ? 1 - %i',size(CONTRAST,1)),1);
	c = CONTRAST(i,:);
end


% get height threshold [default = 3.2]
%---------------------------------------------------------------------------
u    = spm_input('height threshold {Z value}',2,'e',3.2);

% get extent threshold [default = E{n} - expected voxels per cluster]
% Omit spatial extent threshold for multiple contrasts.
%---------------------------------------------------------------------------
if size(c,1) == 1
	[P,EN,Em,En,Pk] = spm_P(1,W,u,0,S);
	k    = spm_input('extent threshold {voxels}',3,'e',round(En));
else
	k    = 0;
end


% pass arguments to spm_projections
%---------------------------------------------------------------------------
set(2,'Name','Thankyou','Pointer','watch')

[t XYZ] = spm_projections(SPMt(i,:),XYZ,u,k,V,W,S,[K H C B G],c,df);


%-Write out filtered SPM{Z} ? (APH addition - 23/06/95)
%-----------------------------------------------------------------------
if size(XYZ,2) %-Only proceed if there's something to work on
    set(2,'Name','SPM{Z} projections','Pointer','arrow')
    if spm_input('Write filtered SPM{Z} to file?',5,'y/n')=='y'
	FName=sprintf('SPM%d_filtered',i);
	FName=spm_input('Filename ?',6,'s',FName);
	str=sprintf('spm{Z}-filtered: corrected p<%6.4f',pV);

	%-Reconstruct filtered image from XYZ & t
	%---------------------------------------------------------------
	n=size(XYZ,2);
	if exist('FLIP')~=1, FLIP=0; end
	%-Unflip flipped images - Negate X values if FLIP==1
	rcp=round(XYZ./meshgrid([1-2*FLIP;1;1].*V(4:6),1:n)' + ...
		meshgrid(V(7:9),1:n)');
	DimMult=cumprod([1,V(1:2)']);
	OffSets=meshgrid([0,1,1],1:n)';
	e=((rcp-OffSets)'*DimMult')';
	t=t.*(t>0); %-Ensure positivity of z-values
	T=zeros(1,prod(V(1:3)));
	T(e)=t;

	%-Write out to analyze file
	%---------------------------------------------------------------
	spm_hwrite([FName,'.hdr'],V(1:3),V(4:6),1/16,2,0,V(7:9),str);
	fid = fopen([FName,'.img'],'w');
	fwrite(fid,T*16,spm_type(2));
	fclose(fid);
    end % (if)
end % (if)

%-Finished
%---------------------------------------------------------------------------
figure(2); clf
set(2,'Name','','Pointer','arrow')



















