% fif2spm(fifName[, spmName[, overwrite]])
%
% fif2spm takes the fiff file specified in fifName and writes it as an SPM
% (.img) file in the file specified as spmName. If spmName is not specified,
% fifName is used, with the extension changed. If overwrite is set to true,
% files are overwritten without warning; if false, files are not (defaults to
% true).
%
% Based on the script fif2spm_script.m by Su Li (?).
%
% CW 9-2010

function fif2spm(fifName, spmName, overwrite)

	if ~exist('spmName', 'var')
		spmName = [fifName(1:end-3) 'img'];
	end%if
	if ~exist('overwrite', 'var')
		overwrite = true;
	end%if
	
	[fifPath fifName fifExt fifVsn] = fileparts(fifName);
	fifName = [fifName fifExt fifVsn];
	clear fifExt fifVsn;
	
	[spmPath spmName spmExt spmVsn] = fileparts(spmName);
	spmName = [spmName spmExt spmVsn];
	clear spmExt spmVsn;	
	
	%%...?

end%function