function spm_append(MAT,X)
% appends a matrix to a MAT-file (Level 1 format)
% FORMAT spm_append(MAT,X);
% MAT  - name of .mat file
% X    - matrix
%____________________________________________________________________________
%
% spm_append is used to append variables to a matrix in a MAT-
% file without loading the matrix into working memory.  This saves
% memory and time.  spm_append is used by spm_spm as the latter
% cycles over planes saving data for subsequent analysis
%
% The name of the matrix is set to MAT ensuring only one matrix per
% MAT-file.  MAT.mat is created if necessary.
%
% NOTE: this routine will require updating for Level 2 MAT-file format and
% forwards compatibility with new versions of MATLAB
% spm_append will only append to MAT-files in pwd
%
%__________________________________________________________________________
% %W% %E%


%----------------------------------------------------------------------------
if ~length(X); return; end
X     = real(X);

% create if necessary
%----------------------------------------------------------------------------
fid   = fopen([pwd '/' MAT '.mat'],'r');
if fid < 0
        eval([MAT ' = X(:,1);']);
        eval(['save ' MAT '.mat ' MAT]);
        spm_append(MAT,X(:,[2:size(X,2)]));
        return
end

% check for number of rows compatibility
%----------------------------------------------------------------------------
[m n] = size(X);
HDR   = fread(fid,5,'int32');
if m ~= HDR(2);  error('  Incompatible sizes');  end
fclose(fid);

% append and update nunmber of columns
%----------------------------------------------------------------------------
fid   = fopen([MAT '.mat'],'a');
fwrite(fid,X(:),'double');
fclose(fid);
fid   = fopen([MAT '.mat'],'r+');
fseek(fid,8,'bof');
fwrite(fid,(HDR(3) + n),'int32');
fclose(fid);
